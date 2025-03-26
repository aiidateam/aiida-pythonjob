import importlib
import inspect
from typing import Any, Callable, Dict, List, Optional, Tuple, Union, _SpecialForm, get_type_hints

from aiida.common.exceptions import NotExistent
from aiida.engine import ExitCode
from aiida.orm import Computer, InstalledCode, User, load_code, load_computer


def get_required_imports(func: Callable) -> Dict[str, set]:
    """Retrieve type hints and the corresponding modules."""
    type_hints = get_type_hints(func)
    imports = {}

    def add_imports(type_hint):
        if isinstance(type_hint, _SpecialForm):  # Handle special forms like Any, Union, Optional
            module_name = "typing"
            type_name = type_hint._name or str(type_hint)
        elif hasattr(type_hint, "__origin__"):  # This checks for higher-order types like List, Dict
            module_name = type_hint.__module__
            type_name = getattr(type_hint, "_name", None) or getattr(type_hint.__origin__, "__name__", None)
            for arg in getattr(type_hint, "__args__", []):
                if arg is type(None):
                    continue
                add_imports(arg)  # Recursively add imports for each argument
        elif hasattr(type_hint, "__module__"):
            module_name = type_hint.__module__
            type_name = type_hint.__name__
        else:
            return  # If no module or origin, we can't import it, e.g., for literals
        if type_name is not None:
            if module_name not in imports:
                imports[module_name] = set()
            imports[module_name].add(type_name)

    for _, type_hint in type_hints.items():
        add_imports(type_hint)
    return imports


def inspect_function(
    func: Callable, inspect_source: bool = False, register_pickle_by_value: bool = False
) -> Dict[str, Any]:
    """Serialize a function for storage or transmission."""
    # we need save the source code explicitly, because in the case of jupyter notebook,
    # the source code is not saved in the pickle file
    import cloudpickle

    if inspect_source:
        try:
            source_code = inspect.getsource(func)
            # Split the source into lines for processing
            source_code_lines = source_code.split("\n")
            source_code = "\n".join(source_code_lines)
        except OSError:
            source_code = "Failed to retrieve source code."
    else:
        source_code = ""

    if register_pickle_by_value:
        module = importlib.import_module(func.__module__)
        cloudpickle.register_pickle_by_value(module)
        pickled_function = cloudpickle.dumps(func)
        cloudpickle.unregister_pickle_by_value(module)
    else:
        pickled_function = cloudpickle.dumps(func)

    return {"source_code": source_code, "mode": "use_pickled_function", "pickled_function": pickled_function}


def build_function_data(func: Callable, register_pickle_by_value: bool = False) -> Dict[str, Any]:
    """Inspect the function and return a dictionary with the function data."""
    import types

    if isinstance(func, (types.FunctionType, types.BuiltinFunctionType, type)):
        # Check if callable is nested (contains dots in __qualname__ after the first segment)
        function_data = {"name": func.__name__}
        if func.__module__ == "__main__" or "." in func.__qualname__.split(".", 1)[-1]:
            # Local or nested callable, so pickle the callable
            function_data.update(inspect_function(func, inspect_source=True))
        else:
            # Global callable (function/class), store its module and name for reference
            function_data.update(inspect_function(func, register_pickle_by_value=register_pickle_by_value))
    else:
        raise TypeError("Provided object is not a callable function or class.")
    return function_data


def get_or_create_code(
    label: str = "python3",
    computer: Optional[Union[str, "Computer"]] = "localhost",
    filepath_executable: Optional[str] = None,
    prepend_text: str = "",
) -> InstalledCode:
    """Try to load code, create if not exit."""

    try:
        return load_code(f"{label}@{computer}")
    except NotExistent:
        description = f"Code on computer: {computer}"
        computer = load_computer(computer)
        filepath_executable = filepath_executable or label
        code = InstalledCode(
            computer=computer,
            label=label,
            description=description,
            filepath_executable=filepath_executable,
            default_calc_job_plugin="pythonjob.pythonjob",
            prepend_text=prepend_text,
        )

        code.store()
        return code


def generate_bash_to_create_python_env(
    name: str,
    pip: Optional[List[str]] = None,
    conda: Optional[Dict[str, list]] = None,
    modules: Optional[List[str]] = None,
    python_version: Optional[str] = None,
    variables: Optional[Dict[str, str]] = None,
    shell: str = "posix",
):
    """
    Generates a bash script for creating or updating a Python environment on a remote computer.
    If python_version is None, it uses the Python version from the local environment.
    Conda is a dictionary that can include 'channels' and 'dependencies'.
    """
    import sys

    pip = pip or []
    conda_channels = conda.get("channels", []) if conda else []
    conda_dependencies = conda.get("dependencies", []) if conda else []
    # Determine the Python version from the local environment if not provided
    local_python_version = f"{sys.version_info.major}.{sys.version_info.minor}"
    desired_python_version = python_version if python_version is not None else local_python_version

    # Start of the script
    script = "#!/bin/bash\n\n"

    # Load modules if provided
    if modules:
        script += "# Load specified system modules\n"
        for module in modules:
            script += f"module load {module}\n"

    # Conda shell hook initialization for proper conda activation
    script += "# Initialize Conda for this shell\n"
    script += f'eval "$(conda shell.{shell} hook)"\n'

    script += "# Setup the Python environment\n"
    script += "if ! conda info --envs | grep -q ^{name}$; then\n"
    script += "    # Environment does not exist, create it\n"
    if conda_dependencies:
        dependencies_string = " ".join(conda_dependencies)
        script += f"    conda create -y -n {name} python={desired_python_version} {dependencies_string}\n"
    else:
        script += f"    conda create -y -n {name} python={desired_python_version}\n"
    script += "fi\n"
    if conda_channels:
        script += "EXISTING_CHANNELS=$(conda config --show channels)\n"
        script += "for CHANNEL in " + " ".join(conda_channels) + ";\n"
        script += "do\n"
        script += '    if ! echo "$EXISTING_CHANNELS" | grep -q $CHANNEL; then\n'
        script += "        conda config --prepend channels $CHANNEL\n"
        script += "    fi\n"
        script += "done\n"
    script += f"conda activate {name}\n"

    # Install pip packages
    if pip:
        script += f"pip install {' '.join(pip)}\n"

    # Set environment variables
    if variables:
        for var, value in variables.items():
            script += f"export {var}='{value}'\n"

    # End of the script
    script += "echo 'Environment setup is complete.'\n"

    return script


def create_conda_env(
    computer: Union[str, Computer],
    name: str,
    pip: Optional[List[str]] = None,
    conda: Optional[List[str]] = None,
    modules: Optional[List[str]] = None,
    python_version: Optional[str] = None,
    variables: Optional[Dict[str, str]] = None,
    shell: str = "posix",
) -> Tuple[bool, str]:
    """Test that there is no unexpected output from the connection."""
    # Execute a command that should not return any error, except ``NotImplementedError``
    # since not all transport plugins implement remote command execution.
    from aiida.common.exceptions import NotExistent

    user = User.collection.get_default()
    if isinstance(computer, str):
        computer = load_computer(computer)
    try:
        authinfo = computer.get_authinfo(user)
    except NotExistent:
        raise f"Computer<{computer.label}> is not yet configured for user<{user.email}>"

    scheduler = authinfo.computer.get_scheduler()
    transport = authinfo.get_transport()

    script = generate_bash_to_create_python_env(name, pip, conda, modules, python_version, variables, shell)
    with transport:
        scheduler.set_transport(transport)
        try:
            retval, stdout, stderr = transport.exec_command_wait(script)
        except NotImplementedError:
            return (
                True,
                f"Skipped, remote command execution is not implemented for the "
                f"`{computer.transport_type}` transport plugin",
            )

        if retval != 0:
            return (
                False,
                f"The command returned a non-zero return code ({retval})",
            )

        template = """
We detected an error while creating the environemnt on the remote computer, as shown between the bars
=============================================================================================
{}
=============================================================================================
Please check!
    """
        if stderr:
            return False, template.format(stderr)

        if stdout:
            # the last line is the echo 'Environment setup is complete.'
            if not stdout.strip().endswith("Environment setup is complete."):
                return False, template.format(stdout)
            else:
                return True, "Environment setup is complete."

    return True, None


def format_input_output_ports(data):
    ports = data.get("ports", [])
    if ports:
        data["identifier"] = "NAMESPACE"
        new_ports = []
        for item in ports:
            if isinstance(item, str):
                new_ports.append({"name": item, "identifier": "ANY"})
            elif isinstance(item, dict):
                item.setdefault("identifier", "any")
                # if the output is WORKGRAPH.NAMESPACE, we need to change it to NAMESPACE
                if item["identifier"].split(".")[-1].upper() == "NAMESPACE":
                    item["identifier"] = "NAMESPACE"
                    new_ports.append(format_input_output_ports(item))
                else:
                    new_ports.append(item)
            else:
                raise ValueError(f"Invalid schema: {item}")
        data["ports"] = new_ports
    else:
        data.setdefault("identifier", "ANY")
    return data


def build_input_port_definitions(
    func,
    input_ports: Optional[List[Dict[str, Any]]] = None,
) -> List[Dict[str, Any]]:
    """
    Build a list of port definitions from a function signature, then merge
    user-defined ports that may override or add to them.

    Return: A list of dicts, each at least has:
      - "name":  str
      - "identifier": str   (e.g. "ANY" or "NAMESPACE")
      - possibly other keys like "required", etc.
    """
    import inspect

    signature = inspect.signature(func)
    default_ports = []

    for param_name, param in signature.parameters.items():
        if param.kind in (
            inspect.Parameter.POSITIONAL_ONLY,
            inspect.Parameter.POSITIONAL_OR_KEYWORD,
            inspect.Parameter.KEYWORD_ONLY,
        ):
            default_ports.append(
                {"name": param_name, "identifier": "ANY", "required": param.default is inspect.Parameter.empty}
            )
        elif param.kind == inspect.Parameter.VAR_POSITIONAL:
            # This is *args
            raise NotImplementedError("Variable positional arguments are not yet supported")
        elif param.kind == inspect.Parameter.VAR_KEYWORD:
            # This is **kwargs
            default_ports.append({"name": param_name, "identifier": "NAMESPACE", "required": False})
        else:
            raise ValueError(f"Unsupported parameter kind: {param.kind}")

    # Now merge in user-defined overrides or additions
    user_defined_ports = input_ports.get("ports", []) if input_ports else []

    # Convert default list to a dict by name for easy merging
    merged_dict = {port["name"]: port for port in default_ports}
    for user_port in user_defined_ports:
        name = user_port["name"]
        merged_dict[name] = {**merged_dict.get(name, {}), **user_port}

    # Return the final merged list
    input_ports["ports"] = list(merged_dict.values())
    return input_ports


def serialize_ports(
    python_data: Dict[str, Any],
    port_schema: Dict[str, Any],
    serializers: Optional[Dict[str, str]] = None,
) -> Dict[str, Any]:
    """
    Convert a raw Python dictionary of user inputs into a dictionary of AiiDA nodes,
    according to a nested port_schema.

    :param python_data:   The actual user data, e.g. {"my_namespace": {...}, "some_other_port": ...}.
    :param port_schema:   A dict of port definition, may contain sub ports with each is a dict with at least:
                          { "name": <str>, "identifier": "ANY" or "NAMESPACE", ... }
                          If NAMESPACE, it may have "ports": <list of nested port definitions>.
    :param serializers:   (Optional) custom mapping for specialized serialization.

    :return: A dictionary with the same structure, but with AiiDA nodes (or nested dicts)
             instead of raw Python values.
    """
    from aiida_pythonjob.data.serializer import general_serializer

    if port_schema["identifier"].upper() == "NAMESPACE":
        name = port_schema["name"]
        if not isinstance(python_data, dict):
            raise ValueError(f"Expected dict for namespace '{name}', got {type(python_data)}")
        # Convert schema into a dict keyed by port name for convenience.
        sub_port_map = {p["name"]: p for p in port_schema.get("ports", [])}
        result = {}
        for name, value in python_data.items():
            sub_port = sub_port_map.get(name, {})
            port_id = sub_port.get("identifier", "ANY").upper()
            if port_id == "NAMESPACE":
                result[name] = serialize_ports(value, sub_port, serializers=serializers)

            else:
                serialized = general_serializer(value, serializers=serializers, store=False)
                result[name] = serialized
    else:
        serialized = general_serializer(python_data, serializers=serializers, store=False)
        result = serialized

    return result


def deserialize_ports(
    serialized_data: Dict[str, Any],
    port_schema: Dict[str, Any],
    deserializers: Optional[Dict[str, str]] = None,
) -> Dict[str, Any]:
    """
    Convert a dictionary of AiiDA data nodes (or nested dicts) back into raw Python objects,
    according to the same nested port_schema.

    :param serialized_data: The data stored in AiiDA, e.g.
                            { "my_namespace": { "foo": <Int>, "bar": { "baz": <List> } }, ... }
    :param port_schema:     The same schema that was used to serialize.
    :param deserializers:   (Optional) custom mapping for specialized deserialization.

    :return: A dictionary of raw Python values.
    """
    from aiida_pythonjob.data.deserializer import deserialize_to_raw_python_data

    if port_schema["identifier"].upper() == "NAMESPACE":
        name = port_schema["name"]
        if not isinstance(serialized_data, dict):
            raise ValueError(f"Expected dict for namespace '{name}', got {type(serialized_data)}")
        sub_port_map = {p["name"]: p for p in port_schema.get("ports", [])}
        result = {}
        for name, node_or_subdict in serialized_data.items():
            sub_port = sub_port_map.get(name, {})
            port_id = sub_port.get("identifier", "ANY").upper()

            if port_id == "NAMESPACE":
                result[name] = deserialize_ports(node_or_subdict, sub_port, deserializers=deserializers)
            else:
                raw_value = deserialize_to_raw_python_data(node_or_subdict, deserializers=deserializers)
                result[name] = raw_value
    else:
        raw_value = deserialize_to_raw_python_data(serialized_data, deserializers=deserializers)
        result = raw_value

    return result


def already_serialized(results):
    """Check if the results are already serialized."""
    import collections

    from aiida import orm

    if isinstance(results, orm.Data):
        return True
    elif isinstance(results, collections.abc.Mapping):
        for value in results.values():
            if not already_serialized(value):
                return False
        return True
    else:
        return False


def parse_outputs(
    results: Any,
    output_ports: Dict[str, Any],
    exit_codes,
    logger,
    serializers: Optional[Dict[str, str]] = None,
) -> Union[Dict[str, Any], ExitCode]:
    """
    Parse the "results" returned from a function or loaded from file,
    given a schema of 'output_ports'."
    """

    # Read output_ports specification

    if isinstance(results, tuple):
        if len(output_ports["ports"]) != len(results):
            return exit_codes.ERROR_RESULT_OUTPUT_MISMATCH
        for i in range(len(output_ports["ports"])):
            output_ports["ports"][i]["value"] = serialize_ports(
                python_data=results[i], port_schema=output_ports["ports"][i], serializers=serializers
            )
    elif isinstance(results, dict):
        # pop the exit code if it exists inside the dictionary
        exit_code = results.pop("exit_code", None)
        if exit_code:
            # If there's an exit_code, handle it (dict or int)
            if isinstance(exit_code, dict):
                exit_code = ExitCode(exit_code["status"], exit_code["message"])
            elif isinstance(exit_code, int):
                exit_code = ExitCode(exit_code)
            if exit_code.status != 0:
                return exit_code
        if len(output_ports["ports"]) == 1:
            # User returned a single (nested) dict with AiiDA data nodes as values
            if already_serialized(results):
                output_ports["ports"] = [{"name": key, "value": value} for key, value in results.items()]
            elif output_ports["ports"][0]["name"] in results:
                output_ports["ports"][0]["value"] = serialize_ports(
                    python_data=results.pop(output_ports["ports"][0]["name"]),
                    port_schema=output_ports["ports"][0],
                    serializers=serializers,
                )
                # If there are any extra keys in `results`, log a warning
                if len(results) > 0:
                    logger.warning(f"Found extra results that are not included in the output: {results.keys()}")
            else:
                # Otherwise assume the entire dict is the single output
                output_ports["ports"][0]["value"] = serialize_ports(
                    python_data=results, port_schema=output_ports["ports"][0], serializers=serializers
                )
        elif len(output_ports["ports"]) > 1:
            # Match each top-level output by name
            for output in output_ports["ports"]:
                if output["name"] not in results:
                    if output.get("required", True):
                        return exit_codes.ERROR_MISSING_OUTPUT
                else:
                    output["value"] = serialize_ports(
                        python_data=results.pop(output["name"]), port_schema=output, serializers=serializers
                    )
            # Any remaining results are unaccounted for -> log a warning
            if len(results) > 0:
                logger.warning(f"Found extra results that are not included in the output: {results.keys()}")
    elif len(output_ports["ports"]) == 1:
        # Single top-level output
        # There are two cases:
        # 1. The output as a whole will be serialized as the single output
        # 2. The output is a mapping with already AiiDA data nodes as values, no need to serialize
        if already_serialized(results):
            output_ports["ports"][0]["value"] = results
        else:
            output_ports["ports"][0]["value"] = serialize_ports(
                python_data=results, port_schema=output_ports["ports"][0], serializers=serializers
            )
    else:
        return exit_codes.ERROR_RESULT_OUTPUT_MISMATCH
