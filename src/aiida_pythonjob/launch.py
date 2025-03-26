from __future__ import annotations

import inspect
import os
from typing import Any, Callable, Dict, List, Optional, Union

from aiida import orm

from .utils import (
    build_function_data,
    build_input_port_definitions,
    format_input_output_ports,
    get_or_create_code,
)


def prepare_pythonjob_inputs(
    function: Optional[Callable[..., Any]] = None,
    function_inputs: Optional[Dict[str, Any]] = None,
    input_ports: Optional[List[str | dict]] = None,
    output_ports: Optional[List[str | dict]] = None,
    code: Optional[orm.AbstractCode] = None,
    command_info: Optional[Dict[str, str]] = None,
    computer: Union[str, orm.Computer] = "localhost",
    metadata: Optional[Dict[str, Any]] = None,
    upload_files: Dict[str, str] = {},
    process_label: Optional[str] = None,
    function_data: dict | None = None,
    deserializers: dict | None = None,
    serializers: dict | None = None,
    register_pickle_by_value: bool = False,
    **kwargs: Any,
) -> Dict[str, Any]:
    """Prepare the inputs for PythonJob"""
    from .utils import serialize_ports

    if function is None and function_data is None:
        raise ValueError("Either function or function_data must be provided")
    if function is not None and function_data is not None:
        raise ValueError("Only one of function or function_data should be provided")
    # if function is a function, inspect it and get the source code
    if function is not None and inspect.isfunction(function):
        function_data = build_function_data(function, register_pickle_by_value=register_pickle_by_value)
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
    output_ports = output_ports or [{"name": "result"}]
    output_ports = {"name": "outputs", "identifier": "namespace", "ports": output_ports or [{"name": "result"}]}
    input_ports = {"name": "inputs", "identifier": "namespace", "ports": input_ports or []}
    input_ports = format_input_output_ports(input_ports)
    input_ports = build_input_port_definitions(func=function, input_ports=input_ports)
    function_data["output_ports"] = format_input_output_ports(output_ports)
    function_data["input_ports"] = input_ports
    # serialize the kwargs into AiiDA Data
    function_inputs = function_inputs or {}
    function_inputs = serialize_ports(python_data=function_inputs, port_schema=input_ports, serializers=serializers)
    # replace "." with "__dot__" in the keys of a dictionary
    if deserializers:
        deserializers = orm.Dict({k.replace(".", "__dot__"): v for k, v in deserializers.items()})
    if serializers:
        serializers = orm.Dict({k.replace(".", "__dot__"): v for k, v in serializers.items()})
    inputs = {
        "function_data": function_data,
        "code": code,
        "function_inputs": function_inputs,
        "upload_files": new_upload_files,
        "metadata": metadata or {},
        "deserializers": deserializers,
        "serializers": serializers,
        **kwargs,
    }
    if process_label:
        inputs["process_label"] = process_label
    return inputs


def create_inputs(func, *args: Any, **kwargs: Any) -> dict[str, Any]:
    """Create the input dictionary for the ``FunctionProcess``."""
    # The complete input dictionary consists of the keyword arguments...
    inputs = dict(kwargs or {})
    arguments = list(args)
    for name, parameter in inspect.signature(func).parameters.items():
        if parameter.kind in [parameter.POSITIONAL_ONLY, parameter.POSITIONAL_OR_KEYWORD]:
            try:
                inputs[name] = arguments.pop(0)
            except IndexError:
                pass
        elif parameter.kind is parameter.VAR_POSITIONAL:
            raise NotImplementedError("Variable positional arguments are not yet supported")

    return inputs


def prepare_pyfunction_inputs(
    function: Optional[Callable[..., Any]] = None,
    function_inputs: Optional[Dict[str, Any]] = None,
    input_ports: Optional[List[str | dict]] = None,
    output_ports: Optional[List[str | dict]] = None,
    metadata: Optional[Dict[str, Any]] = None,
    process_label: Optional[str] = None,
    function_data: dict | None = None,
    deserializers: dict | None = None,
    serializers: dict | None = None,
    register_pickle_by_value: bool = False,
    **kwargs: Any,
) -> Dict[str, Any]:
    """Prepare the inputs for PythonJob"""
    import types

    from .utils import serialize_ports

    if function is None and function_data is None:
        raise ValueError("Either function or function_data must be provided")
    if function is not None and function_data is not None:
        raise ValueError("Only one of function or function_data should be provided")
    # if function is a function, inspect it and get the source code
    if function is not None:
        if inspect.isfunction(function):
            function_data = build_function_data(function, register_pickle_by_value=register_pickle_by_value)
        elif isinstance(function, types.BuiltinFunctionType):
            raise NotImplementedError("Built-in functions are not supported yet")
        else:
            raise ValueError("Invalid function type")
    output_ports = {"name": "outputs", "identifier": "namespace", "ports": output_ports or [{"name": "result"}]}
    input_ports = {"name": "inputs", "identifier": "namespace", "ports": input_ports or []}
    input_ports = format_input_output_ports(input_ports)
    input_ports = build_input_port_definitions(func=function, input_ports=input_ports)
    function_data["output_ports"] = format_input_output_ports(output_ports)
    function_data["input_ports"] = input_ports
    # serialize the kwargs into AiiDA Data
    function_inputs = function_inputs or {}
    function_inputs = serialize_ports(python_data=function_inputs, port_schema=input_ports, serializers=serializers)
    # replace "." with "__dot__" in the keys of a dictionary
    if deserializers:
        deserializers = orm.Dict({k.replace(".", "__dot__"): v for k, v in deserializers.items()})
    if serializers:
        serializers = orm.Dict({k.replace(".", "__dot__"): v for k, v in serializers.items()})
    inputs = {
        "function_data": function_data,
        "function_inputs": function_inputs,
        "metadata": metadata or {},
        "deserializers": deserializers,
        "serializers": serializers,
        **kwargs,
    }
    if process_label:
        inputs["process_label"] = process_label
    return inputs
