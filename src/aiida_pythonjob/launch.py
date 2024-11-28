import inspect
import os
from typing import Any, Callable, Dict, Optional, Union

from aiida.orm import AbstractCode, Computer, FolderData, List, SinglefileData, Str

from .data.pickled_function import PickledFunction
from .data.serializer import serialize_to_aiida_nodes
from .utils import get_or_create_code


def prepare_pythonjob_inputs(
    function: Optional[Callable[..., Any]] = None,
    function_inputs: Optional[Dict[str, Any]] = None,
    function_outputs: Optional[Dict[str, Any]] = None,
    code: Optional[AbstractCode] = None,
    command_info: Optional[Dict[str, str]] = None,
    computer: Union[str, Computer] = "localhost",
    metadata: Optional[Dict[str, Any]] = None,
    upload_files: Dict[str, str] = {},
    process_label: Optional[str] = None,
    pickled_function: Optional[PickledFunction] = None,
    **kwargs: Any,
) -> Dict[str, Any]:
    pass
    """Prepare the inputs for PythonJob"""

    if function is None and pickled_function is None:
        raise ValueError("Either function or pickled_function must be provided")
    if function is not None and pickled_function is not None:
        raise ValueError("Only one of function or pickled_function should be provided")
    # if function is a function, convert it to a PickledFunction
    if function is not None and inspect.isfunction(function):
        executor = PickledFunction.build_callable(function)
    if pickled_function is not None:
        executor = pickled_function
    new_upload_files = {}
    # change the string in the upload files to SingleFileData, or FolderData
    for key, source in upload_files.items():
        # only alphanumeric and underscores are allowed in the key
        # replace all "." with "_dot_"
        new_key = key.replace(".", "_dot_")
        if isinstance(source, str):
            if os.path.isfile(source):
                new_upload_files[new_key] = SinglefileData(file=source)
            elif os.path.isdir(source):
                new_upload_files[new_key] = FolderData(tree=source)
        elif isinstance(source, (SinglefileData, FolderData)):
            new_upload_files[new_key] = source
        else:
            raise ValueError(f"Invalid upload file type: {type(source)}, {source}")
    #
    if code is None:
        command_info = command_info or {}
        code = get_or_create_code(computer=computer, **command_info)
    # get the source code of the function
    function_name = executor["name"]
    if executor.get("is_pickle", False):
        function_source_code = executor["import_statements"] + "\n" + executor["source_code_without_decorator"]
    else:
        function_source_code = f"from {executor['module']} import {function_name}"

    # serialize the kwargs into AiiDA Data
    function_inputs = function_inputs or {}
    function_inputs = serialize_to_aiida_nodes(function_inputs)
    # transfer the args to kwargs
    inputs = {
        "process_label": process_label or "PythonJob<{}>".format(function_name),
        "function_source_code": Str(function_source_code),
        "function_name": Str(function_name),
        "code": code,
        "function_inputs": function_inputs,
        "upload_files": new_upload_files,
        "function_outputs": List(function_outputs),
        "metadata": metadata or {},
        **kwargs,
    }
    return inputs
