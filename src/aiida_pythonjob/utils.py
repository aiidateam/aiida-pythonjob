import inspect
from typing import Any, Callable, Dict, List, Optional, Tuple, Union, _SpecialForm, get_type_hints

from aiida.common.exceptions import NotExistent
from aiida.orm import Computer, InstalledCode, User, load_code, load_computer

CONDA_DEFAULT_PATH = "$HOME/miniforge3/"


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


def inspect_function(func: Callable) -> Dict[str, Any]:
    """Serialize a function for storage or transmission."""
    # we need save the source code explicitly, because in the case of jupyter notebook,
    # the source code is not saved in the pickle file
    from aiida_pythonjob.data.pickled_data import PickledData

    try:
        source_code = inspect.getsource(func)
        # Split the source into lines for processing
        source_code_lines = source_code.split("\n")
        source_code = "\n".join(source_code_lines)
    except OSError:
        source_code = "Failed to retrieve source code."

    return {"source_code": source_code, "mode": "use_pickled_function", "pickled_function": PickledData(value=func)}


def build_function_data(func: Callable) -> Dict[str, Any]:
    """Inspect the function and return a dictionary with the function data."""
    import types

    if isinstance(func, (types.FunctionType, types.BuiltinFunctionType, type)):
        # Check if callable is nested (contains dots in __qualname__ after the first segment)
        function_data = {"name": func.__name__}
        if func.__module__ == "__main__" or "." in func.__qualname__.split(".", 1)[-1]:
            # Local or nested callable, so pickle the callable
            function_data.update(inspect_function(func))
        else:
            # Global callable (function/class), store its module and name for reference
            function_data.update(
                {
                    "mode": "use_module_path",
                    "source_code": f"from {func.__module__} import {func.__name__}",
                }
            )
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


def generate_bash_to_install_conda(
    shell: str = "posix",
    destination: str = CONDA_DEFAULT_PATH,
    modules: Optional[list] = None,
):
    """
    Args:
        shell (str): The type of shell to initialize conda for (default is "posix").
        destination (str): The installation directory for Miniforge (default is CONDA_DEFAULT_PATH).
        modules (list): A list of system modules to load before running the script (default is None).
    Returns:
        str: A bash script as a string to install Miniforge and set up conda.

    Generates a bash script to install conda via miniforge on a local/remote computer.
    The default channel (the only one) is automatically set to be conda-forge, avoiding then to
    use Anaconda channels, restricted by the license.
    We anyway perform a check to be sure that the installation will not use Anaconda channels.
    If python_version is None, it uses the Python version from the local environment.
    """

    # Start of the script
    script = "#!/bin/bash\n\n"

    # Load modules if provided
    if modules:
        script += "# Load specified system modules\n"
        for module in modules:
            script += f"module load {module}\n"

    script += f"""
# Check if conda is already installed
if command -v {destination}/bin/conda &> /dev/null; then
    echo "Conda is already installed. Skipping installation."
else\n
"""

    # Getting minimum Miniforge installer as recommended here: https://github.com/conda-forge/miniforge?tab=readme-ov-file
    script += "# Downloading Miniforge installer\n"
    script += "curl -L -O \
        https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh\n"

    # Running the installer
    script += "# Running the Miniforge installer\n"
    script += f"bash Miniforge3-$(uname)-$(uname -m).sh -b -p {destination}\n"

    # Conda shell hook initialization for proper conda activation
    script += "# Initialize Conda for this shell\n"
    script += f'eval "$({destination}/bin/conda shell.{shell} hook)"\n'

    # Ensure the default Anaconda channel is not present in the conda configuration
    script += """
# Check if 'conda config --show channels | grep default' returns anything
if conda config --show channels | grep -q "defaults"; then
    echo "The default Anaconda channel is present in the conda configuration. We remove it."
    conda config --remove channels defaults
else
    echo "The default Anaconda channel is not present in the conda configuration. Good."
fi
"""

    # Ensure the conda-forge channel is present in the conda configuration
    script += """
# Ensure conda-forge is there
if conda config --show channels | grep -q "conda-forge"; then
    echo "The conda-forge channel is present in the conda configuration. Good."
else
    echo "The conda-forge channel is not present in the conda configuration. We add it."
    conda config --append channels conda-forge
fi
"""

    script += "fi\n"
    # End of the script
    script += 'echo "Miniforge-based conda installation is complete."\n\n'

    return script


def generate_bash_to_create_python_env(
    name: str,
    pip: Optional[List[str]] = None,
    conda: Optional[Dict[str, list]] = {},
    modules: Optional[List[str]] = None,
    python_version: Optional[str] = None,
    variables: Optional[Dict[str, str]] = None,
    shell: str = "posix",
):
    """
    Generates a bash script for creating or updating a Python environment on a remote computer.
    If python_version is None, it uses the Python version from the local environment.
    Conda is a dictionary that can include 'channels' and 'dependencies' and 'path', where 'path' is the path to the
    conda executable (not included in the path), and is needed only to activate the environment.
    """
    import sys

    pip = pip or []
    conda_channels = conda.get("channels", []) if conda else []
    conda_dependencies = conda.get("dependencies", []) if conda else []
    conda_path = conda.get("path", CONDA_DEFAULT_PATH)
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
    script += f'eval "$({conda_path}/bin/conda shell.{shell} hook)"\n'

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
    conda: Optional[Dict[str, list]] = {},
    modules: Optional[List[str]] = None,
    python_version: Optional[str] = None,
    variables: Optional[Dict[str, str]] = None,
    shell: str = "posix",
    install_conda: bool = False,
) -> Tuple[bool, str]:
    """

    Create a conda environment on a remote computer.

    Parameters:
    - computer (Union[str, Computer]): The computer on which to create the environment.
        Can be a string (computer label) or a Computer object.
    - name (str): The name of the conda environment to create.
    - pip (Optional[List[str]]): List of pip packages to install in the environment.
    - conda (Optional[List[str]]): List of conda packages to install in the environment. See the
        `generate_bash_to_create_python_env` function for details.
    - modules (Optional[List[str]]): List of modules to load before creating the environment.
    - python_version (Optional[str]): The Python version to use for the environment.
    - variables (Optional[Dict[str, str]]): Environment variables to set during the environment creation.
    - shell (str): The shell type to use (default is "posix").
    - install_conda (bool): Whether to install conda if it is not already installed (default is False).

    Returns:
    - Tuple[bool, str]: A tuple containing a boolean indicating success or failure, and a string message with details.

    Test that there is no unexpected output from the connection."""
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

    conda_path = conda.get("path", CONDA_DEFAULT_PATH)

    if install_conda:
        install_conda_script = generate_bash_to_install_conda(shell, destination=conda_path, modules=modules)
        script = install_conda_script + script

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
