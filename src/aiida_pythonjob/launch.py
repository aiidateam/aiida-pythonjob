from __future__ import annotations

import inspect
import os
from typing import Any, Callable, Dict, List, Optional, Union

from aiida import orm

from .data.serializer import serialize_to_aiida_nodes
from .utils import build_function_data, get_or_create_code


def prepare_pythonjob_inputs(
    function: Optional[Callable[..., Any]] = None,
    function_inputs: Optional[Dict[str, Any]] = None,
    function_outputs: Optional[List[str | dict]] = None,
    code: Optional[orm.AbstractCode] = None,
    command_info: Optional[Dict[str, str]] = None,
    computer: Union[str, orm.Computer] = "localhost",
    metadata: Optional[Dict[str, Any]] = None,
    upload_files: Dict[str, str] = {},
    process_label: Optional[str] = None,
    function_data: dict | None = None,
    **kwargs: Any,
) -> Dict[str, Any]:
    pass
    """Prepare the inputs for PythonJob"""

    if function is None and function_data is None:
        raise ValueError("Either function or function_data must be provided")
    if function is not None and function_data is not None:
        raise ValueError("Only one of function or function_data should be provided")
    # if function is a function, inspect it and get the source code
    if function is not None and inspect.isfunction(function):
        function_data = build_function_data(function)
    new_upload_files = {}
    # change the string in the upload files to SingleFileData, or FolderData
    for key, source in upload_files.items():
        # only alphanumeric and underscores are allowed in the key
        # replace all "." with "_dot_"
        new_key = key.replace(".", "_dot_")
        if isinstance(source, str):
            if os.path.isfile(source):
                new_upload_files[new_key] = orm.SinglefileData(file=source)
            elif os.path.isdir(source):
                new_upload_files[new_key] = orm.FolderData(tree=source)
            else:
                raise ValueError(f"Invalid upload file path: {source}")
        elif isinstance(source, (orm.SinglefileData, orm.FolderData)):
            new_upload_files[new_key] = source
        else:
            raise ValueError(f"Invalid upload file type: {type(source)}, {source}")
    #
    if code is None:
        command_info = command_info or {}
        code = get_or_create_code(computer=computer, **command_info)
    # get the source code of the function
    function_name = function_data["name"]
    if function_data.get("is_pickle", False):
        function_source_code = (
            function_data["import_statements"] + "\n" + function_data["source_code_without_decorator"]
        )
    else:
        function_source_code = f"from {function_data['module']} import {function_name}"

    # serialize the kwargs into AiiDA Data
    function_inputs = function_inputs or {}
    function_inputs = serialize_to_aiida_nodes(function_inputs)
    # transfer the args to kwargs
    inputs = {
        "process_label": process_label or "PythonJob<{}>".format(function_name),
        "function_data": orm.Dict(
            {
                "source_code": function_source_code,
                "name": function_name,
                "outputs": function_outputs or [],
            }
        ),
        "code": code,
        "function_inputs": function_inputs,
        "upload_files": new_upload_files,
        "metadata": metadata or {},
        **kwargs,
    }
    return inputs
