import importlib
import inspect
from typing import (
    Any,
    Callable,
    Dict,
    List,
    Optional,
    Tuple,
    Union,
    _SpecialForm,
    get_type_hints,
)

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


def serialize_ports(
    python_data: Any,
    port_schema: Dict[str, Any],
    serializers: Optional[Dict[str, str]] = None,
) -> Any:
    from aiida_pythonjob.data.serializer import general_serializer

    if port_schema["identifier"].upper() != "NAMESPACE":
        return general_serializer(python_data, serializers=serializers, store=False)

    # Namespace
    name = port_schema.get("name", "<namespace>")
    if not isinstance(python_data, dict):
        raise ValueError(f"Expected dict for namespace '{name}', got {type(python_data)}")

    sub_ports: Dict[str, Dict[str, Any]] = port_schema.get("ports", {}) or {}
    is_dyn = bool(port_schema.get("dynamic"))
    item_schema = port_schema.get("item") if is_dyn else None

    out: Dict[str, Any] = {}
    for key, value in python_data.items():
        if key in sub_ports:
            sub = sub_ports[key]
            if sub.get("identifier", "ANY").upper() == "NAMESPACE":
                out[key] = serialize_ports(value, sub, serializers=serializers)
            else:
                out[key] = general_serializer(value, serializers=serializers, store=False)
        elif is_dyn and item_schema is not None:
            if item_schema.get("identifier", "ANY").upper() == "NAMESPACE":
                out[key] = serialize_ports(value, item_schema, serializers=serializers)
            else:
                out[key] = general_serializer(value, serializers=serializers, store=False)
        else:
            raise ValueError(f"Unexpected key '{key}' for namespace '{name}' (not dynamic).")

    return out


def deserialize_ports(
    serialized_data: Any,
    port_schema: Dict[str, Any],
    deserializers: Optional[Dict[str, str]] = None,
) -> Any:
    from aiida_pythonjob.data.deserializer import deserialize_to_raw_python_data

    if port_schema["identifier"].upper() != "NAMESPACE":
        return deserialize_to_raw_python_data(serialized_data, deserializers=deserializers)

    # Namespace
    name = port_schema.get("name", "<namespace>")
    if not isinstance(serialized_data, dict):
        raise ValueError(f"Expected dict for namespace '{name}', got {type(serialized_data)}")

    sub_ports: Dict[str, Dict[str, Any]] = port_schema.get("ports", {}) or {}
    is_dyn = bool(port_schema.get("dynamic"))
    item_schema = port_schema.get("item") if is_dyn else None

    out: Dict[str, Any] = {}
    for key, value in serialized_data.items():
        if key in sub_ports:
            sub = sub_ports[key]
            if sub.get("identifier", "ANY").upper() == "NAMESPACE":
                out[key] = deserialize_ports(value, sub, deserializers=deserializers)
            else:
                out[key] = deserialize_to_raw_python_data(value, deserializers=deserializers)
        elif is_dyn and item_schema is not None:
            if item_schema.get("identifier", "ANY").upper() == "NAMESPACE":
                out[key] = deserialize_ports(value, item_schema, deserializers=deserializers)
            else:
                out[key] = deserialize_to_raw_python_data(value, deserializers=deserializers)
        else:
            raise ValueError(f"Unexpected key '{key}' for namespace '{name}' (not dynamic).")

    return out


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


def _ordered_port_names(ports: Dict[str, Dict[str, Any]]) -> List[str]:
    return list(ports.keys())  # insertion order (Py3.7+)


def parse_outputs(
    results: Any,
    output_ports: Dict[str, Any],
    exit_codes,
    logger,
    serializers: Optional[Dict[str, str]] = None,
) -> Union[Dict[str, Any], ExitCode]:
    """Populate output_ports['ports'][name]['value'] from results based on schema."""
    ports: Dict[str, Dict[str, Any]] = output_ports.get("ports", {}) or {}

    is_dyn = bool(output_ports.get("dynamic"))
    item_schema = output_ports.get("item") if is_dyn else None

    # tuple -> map by order of port names
    if isinstance(results, tuple):
        names = _ordered_port_names(ports)
        if len(names) != len(results):
            return exit_codes.ERROR_RESULT_OUTPUT_MISMATCH
        for i, name in enumerate(names):
            port = ports[name]
            port["value"] = serialize_ports(results[i], port, serializers=serializers)
        return None

    # dict
    if isinstance(results, dict):
        # handle optional inline exit code
        exit_code = results.pop("exit_code", None)
        if exit_code:
            if isinstance(exit_code, dict):
                exit_code = ExitCode(exit_code["status"], exit_code["message"])
            elif isinstance(exit_code, int):
                exit_code = ExitCode(exit_code)
            if exit_code.status != 0:
                return exit_code

        if len(ports) == 1 and not is_dyn:
            # single output:
            # - if user used the same key as port name, use that value;
            # - else treat the entire dict as the value for that single port.
            ((only_name, only_port),) = ports.items()
            if already_serialized(results):
                output_ports["ports"] = {key: {"value": results[key]} for key in results}
            elif only_name in results:
                only_port["value"] = serialize_ports(results.pop(only_name), only_port, serializers=serializers)
                if results:
                    logger.warning(f"Found extra results that are not included in the output: {list(results.keys())}")
            else:
                only_port["value"] = serialize_ports(results, only_port, serializers=serializers)
            return None

        # multi output:
        #  - match fixed outputs by name
        #  - if top-level is dynamic, treat remaining keys as dynamic items using item_schema
        remaining = dict(results)
        # fixed fields
        for name, port in ports.items():
            if name in remaining:
                port["value"] = serialize_ports(remaining.pop(name), port, serializers=serializers)
            elif port.get("required", True):
                logger.warning(f"Missing required output: {name}")
                return exit_codes.ERROR_MISSING_OUTPUT
        # dynamic items at top-level
        if is_dyn:
            if item_schema is None:
                logger.warning("Outputs marked dynamic but missing 'item' schema; treating as ANY.")
                item_schema = {"identifier": "ANY"}
            for name, value in remaining.items():
                # Create a value entry for the dynamic key directly at top-level
                output_ports["ports"][name] = {"value": serialize_ports(value, item_schema, serializers=serializers)}
            return None
        # not dynamic → leftovers are unexpected
        if remaining:
            logger.warning(f"Found extra results that are not included in the output: {list(remaining.keys())}")
        return None

    # single output + non-dict/tuple
    if len(ports) == 1:
        # Single top-level output
        # There are two cases:
        # 1. The output as a whole will be serialized as the single output
        # 2. The output is a mapping with already AiiDA data nodes as values, no need to serialize
        ((only_name, only_port),) = ports.items()
        if already_serialized(results):
            output_ports["ports"] = {key: {"value": results[key]} for key in results}
        else:
            only_port["value"] = serialize_ports(results, only_port, serializers=serializers)
            return None

    return exit_codes.ERROR_RESULT_OUTPUT_MISMATCH
