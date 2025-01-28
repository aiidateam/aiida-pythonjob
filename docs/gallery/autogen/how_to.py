"""
How to guides
===============

"""


######################################################################
# Preparing inputs for `PythonJob`
# --------------------------------
# The `prepare_pythonjob_inputs` function is available for setting up the
# inputs for a `PythonJob` calculation. This function simplifies the process
# of preparing and serializing data, and configuring the execution environment.
#
# - **Code**: You can specify the `computer` where the job will run, which will
#   create a `python3@computer` code if it doesn't already exist. Alternatively,
#   if the code has already been created, you can set the `code` directly.
#
# - **Data**: Use standard Python data types for input. The `prepare_pythonjob_inputs`
#   function handles the conversion to AiiDA data. For serialization:
#    - The function first searches for an AiiDA data entry point corresponding to the module
#      and class names (e.g., `ase.atoms.Atoms`).
#    - If a matching entry point exists, it is used for serialization.
#    - If no match is found, the data is serialized into binary format using `PickledData`.
#
# - **Python Version**: Ensure the Python version on the remote computer matches the local environment.
#   This is important since pickle is used for data storage and retrieval. Use **conda** to
#   create and activate a virtual environment with the same Python version. Pass metadata
#   to the scheduler to activate the environment during the job execution:
#
#   .. code-block:: python
#
#       metadata = {
#           "options": {
#               'custom_scheduler_commands': 'module load anaconda\nconda activate py3.11\n',
#           }
#       }
#
# --------------------------------------------------
# Create a conda environment on the remote computer
# --------------------------------------------------
# One can use the `create_conda_env` function to create a conda environment
# on the remote computer. The function will create a conda environment with
# the specified packages and modules. The function will update the packages
# if the environment already exists.
#
# .. code-block:: python
#
#     from aiida_pythonjob.utils import create_conda_env
#     # create a conda environment on remote computer
#     create_conda_env(
#            "merlin6",                # Remote computer, already stored in the AiiDA database
#            "test_pythonjob",         # Name of the conda environment
#            modules=["anaconda"],     # Modules to load (e.g., Anaconda)
#            pip=["numpy", "matplotlib"],  # Python packages to install via pip
#            conda={                   # Conda-specific settings
#                "channels": ["conda-forge"],  # Channels to use
#                "dependencies": ["qe"]       # Conda packages to install
#            }
#        )
#
#
# If you don't have conda installed on the remote computer, or you can use Anaconda for license reasons,
# you can run the `create_conda_env` function
# with the `install_conda` parameter set to `True`. This will install conda on the remote computer via the
# miniforge installer. You can find more information about Miniforge and download the installer from the
# official [Miniforge GitHub repository](https://github.com/conda-forge/miniforge). The `conda` dictionary
# can be used to specify the desired conda environment path, adding a new key "path": "/path/to/conda"`.,
# e.g.:
#
# .. code-block:: python
#
#     from aiida_pythonjob.utils import create_conda_env
#     # create a conda environment on remote computer
#     create_conda_env(
#            "merlin6",                # Remote computer, already stored in the AiiDA database
#            "test_pythonjob",         # Name of the conda environment
#            modules=["anaconda"],     # Modules to load (e.g., Anaconda)
#            pip=["numpy", "matplotlib"],  # Python packages to install via pip
#            conda={                   # Conda-specific settings
#                "channels": ["conda-forge"],  # Channels to use
#                "dependencies": ["qe"]       # Conda packages to install
#                "path": "$HOME/miniforge3/" # path to the (new) conda installation
#           },
#        )
#
#
# By default, the conda path will be set to `$HOME/miniforge3/`, if the `path` key is not provided.
#

######################################################################
# Default outputs
# --------------
#
# The default output of the function is `result`. The `PythonJob` task
# will store the result as one node in the database with the key `result`.
#
from aiida import load_profile
from aiida.engine import run_get_node
from aiida_pythonjob import PythonJob, prepare_pythonjob_inputs

load_profile()


def add(x, y):
    return x + y


inputs = prepare_pythonjob_inputs(
    add,
    function_inputs={"x": 1, "y": 2},
    computer="localhost",
)
result, node = run_get_node(PythonJob, inputs=inputs)
print("result: ", result["result"])

######################################################################
# Custom outputs
# --------------
# If the function return a dictionary with fixed number of keys, and you
# want to store the values as separate outputs, you can specify the `function_outputs` parameter.
# For a dynamic number of outputs, you can use the namespace output, which is explained later.
#


def add(x, y):
    return {"sum": x + y, "diff": x - y}


inputs = prepare_pythonjob_inputs(
    add,
    function_inputs={"x": 1, "y": 2},
    function_outputs=[
        {"name": "sum"},
        {"name": "diff"},
    ],
)
result, node = run_get_node(PythonJob, **inputs)

print("result: ")
print("sum: ", result["sum"])
print("diff: ", result["diff"])


######################################################################
# Using parent folder
# --------------
# The parent_folder parameter allows a task to access the output files of
# a parent task. This feature is particularly useful when you want to reuse
# data generated by a previous computation in subsequent computations. In
# the following example, the multiply task uses the `result.txt` file created by the add task.
#
#


def add(x, y):
    z = x + y
    with open("result.txt", "w") as f:
        f.write(str(z))
    return x + y


def multiply(x, y):
    with open("parent_folder/result.txt", "r") as f:
        z = int(f.read())
    return x * y + z


inputs1 = prepare_pythonjob_inputs(
    add,
    function_inputs={"x": 1, "y": 2},
    function_outputs=[{"name": "sum"}],
)

result1, node1 = run_get_node(PythonJob, inputs=inputs1)

inputs2 = prepare_pythonjob_inputs(
    multiply,
    function_inputs={"x": 1, "y": 2},
    function_outputs=[{"name": "product"}],
    parent_folder=result1["remote_folder"],
)

result2, node2 = run_get_node(PythonJob, inputs=inputs2)

print("result: ", result2)

######################################################################
# Upload files or folders to the remote computer
# --------------
# The `upload_files` parameter allows users to upload files or folders to
# the remote computer. The files will be uploaded to the working directory of the remote computer.
#

import os  # noqa: E402

# create a temporary file "input.txt" in the current directory
with open("/tmp/input.txt", "w") as f:
    f.write("2")

# create a temporary folder "inputs_folder" in the current directory
# and add a file "another_input.txt" in the folder
os.makedirs("/tmp/inputs_folder", exist_ok=True)
with open("/tmp/inputs_folder/another_input.txt", "w") as f:
    f.write("3")


def add():
    with open("input.txt", "r") as f:
        a = int(f.read())
    with open("inputs_folder/another_input.txt", "r") as f:
        b = int(f.read())
    return a + b


# ------------------------- Submit the calculation -------------------
# we need use full path to the file
input_file = os.path.abspath("/tmp/input.txt")
input_folder = os.path.abspath("/tmp/inputs_folder")
inputs = prepare_pythonjob_inputs(
    add,
    upload_files={
        "input.txt": input_file,
        "inputs_folder": input_folder,
    },
)
result, node = run_get_node(PythonJob, inputs=inputs)
print("result: ", result["result"])

######################################################################
# Retrieve additional files from the remote computer
# --------------
# Sometimes, one may want to retrieve additional files from the remote
# computer after the job has finished. For example, one may want to retrieve
# the output files generated by the `pw.x` calculation in Quantum ESPRESSO.
#
# One can use the `additional_retrieve_list` parameter to specify which files
# should be retrieved from the working directory and stored in the local
# repository after the job has finished
#


def add(x, y):
    z = x + y
    with open("result.txt", "w") as f:
        f.write(str(z))
    return x + y


inputs = prepare_pythonjob_inputs(
    add,
    function_inputs={"x": 1, "y": 2},
    metadata={
        "options": {
            "additional_retrieve_list": ["result.txt"],
        }
    },
)

result, node = run_get_node(PythonJob, inputs=inputs)
print("retrieved files: ", result["retrieved"].list_object_names())

######################################################################
# Namespace Output
# --------------
#
# The `PythonJob` allows users to define namespace outputs. A namespace output
# is a dictionary with keys and values returned by a function. Each value in
# this dictionary will be serialized to AiiDA data, and the key-value pair
# will be stored in the database.
# Why Use Namespace Outputs?
#
# - **Dynamic and Flexible**: The keys and values in the namespace output are
# not fixed and can change based on the task's execution.
# - **Querying**: The data in the namespace output is stored as an AiiDA data
# node, allowing for easy querying and retrieval.
# - **Data Provenance**: When the data is used as input for subsequent tasks,
# the origin of data is tracked.
#
# For example: Consider a molecule adsorption calculation where the namespace
# output stores the surface slabs of the molecule adsorbed on different surface
# sites. The number of surface slabs can vary depending on the surface. These
# output surface slabs can be utilized as input to the next task to calculate the energy.

from ase import Atoms  # noqa: E402
from ase.build import bulk  # noqa: E402


def generate_structures(structure: Atoms, factor_lst: list) -> dict:
    """Scale the structure by the given factor_lst."""
    scaled_structures = {}
    for i in range(len(factor_lst)):
        atoms = structure.copy()
        atoms.set_cell(atoms.cell * factor_lst[i], scale_atoms=True)
        scaled_structures[f"s_{i}"] = atoms
    return {"scaled_structures": scaled_structures}


inputs = prepare_pythonjob_inputs(
    generate_structures,
    function_inputs={"structure": bulk("Al"), "factor_lst": [0.95, 1.0, 1.05]},
    function_outputs=[{"name": "scaled_structures", "identifier": "namespace"}],
)

result, node = run_get_node(PythonJob, inputs=inputs)
print("scaled_structures: ")
for key, value in result["scaled_structures"].items():
    print(key, value)


######################################################################
# What if my calculation fails?
# --------------------------------
#
# The `PythonJobParser` can return specialized exit codes when different
# kinds of errors occur during the calculation:
#
# - ``ERROR_READING_OUTPUT_FILE`` (310):
#     The retrieved output file (e.g., `results.pickle`) could not be opened or read.
#
# - ``ERROR_INVALID_OUTPUT`` (320):
#     The output file is corrupt or contains unexpected/invalid data structures.
#
# - ``ERROR_RESULT_OUTPUT_MISMATCH`` (321):
#     The number of actual results does not match the number/structure of expected outputs.
#
# - ``ERROR_IMPORT_CLOUDPICKLE_FAILED`` (322):
#     The script on the remote machine failed to import `cloudpickle`.
#     The script writes ``error.json`` describing the ImportError.
#
# - ``ERROR_UNPICKLE_INPUTS_FAILED`` (323):
#     The script failed to unpickle the input data (e.g., `inputs.pickle`).
#
# - ``ERROR_UNPICKLE_FUNCTION_FAILED`` (324):
#     The script failed to load the pickled function (e.g., `function.pkl`).
#
# - ``ERROR_FUNCTION_EXECUTION_FAILED`` (325):
#     An exception was raised during the function call.
#
# - ``ERROR_PICKLE_RESULTS_FAILED`` (326):
#     The script failed to pickle (serialize) the final results.
#
# - ``ERROR_SCRIPT_FAILED`` (327):
#     A catch-all exit code if none of the above match. Indicates an unknown/unrecognized error.
#


######################################################################
# Exit Code
# --------------
# Users can define custom exit codes to indicate the status of the task.
#
# When the function returns a dictionary with an `exit_code` key, the system
# automatically parses and uses this code to indicate the task's status. In
# the case of an error, the non-zero `exit_code` value helps identify the specific problem.
#
#


def add(x, y):
    sum = x + y
    if (sum < 0).any():
        exit_code = {"status": 410, "message": "Some elements are negative"}
        return {"sum": sum, "exit_code": exit_code}
    return {"sum": sum}


inputs = prepare_pythonjob_inputs(
    add,
    function_inputs={"x": 1, "y": -21},
)

result, node = run_get_node(PythonJob, inputs=inputs)
print("exit_status:", node.exit_status)
print("exit_message:", node.exit_message)


######################################################################
# Define your data serializer and deserializer
# --------------
#
# PythonJob search data serializer from the `aiida.data` entry point by the
# module name and class name (e.g., `ase.atoms.Atoms`).
#
# In order to let the PythonJob find the serializer, you must register the
# AiiDA data with the following format:
#
# .. code-block:: ini
#
#    [project.entry-points."aiida.data"]
#    abc.ase.atoms.Atoms = "abc.xyz:MyAtomsData"
#
# This will register a data serializer for `ase.atoms.Atoms` data. `abc` is
# the plugin name, the module name is `xyz`, and the AiiDA data class name is
# `AtomsData`. Learn how to create an AiiDA data class `here <https://aiida.readthedocs.io/projects/aiida-core/en/stable/topics/data_types.html#adding-support-for-custom-data-types>`_.
#
# *Avoid duplicate data serializer*: If you have multiple plugins that
# register the same data serializer, the PythonJob will raise an error.
# You can avoid this by selecting the plugin that you want to use in the configuration file.
#
#
# .. code-block:: json
#
#    {
#        "serializers": {
#            "ase.atoms.Atoms": "abc.ase.atoms.AtomsData" # use the full path to the serializer
#        }
#    }
#
# Save the configuration file as `pythonjob.json` in the aiida configuration
# directory (by default, `~/.aiida` directory).
#
# If you want to pass AiiDA Data node as input, and the node does not have a `value` attribute,
# then one must provide a deserializer for it.
#

from aiida import orm  # noqa: E402


def make_supercell(structure, n=2):
    return structure * [n, n, n]


structure = orm.StructureData(cell=[[1, 0, 0], [0, 1, 0], [0, 0, 1]])
structure.append_atom(position=(0.0, 0.0, 0.0), symbols="Li")

inputs = prepare_pythonjob_inputs(
    make_supercell,
    function_inputs={"structure": structure},
    deserializers={
        "aiida.orm.nodes.data.structure.StructureData": "aiida_pythonjob.data.deserializer.structure_data_to_atoms"
    },
)
result, node = run_get_node(PythonJob, inputs=inputs)
print("result: ", result["result"])

######################################################################
# One can also set the deserializer in the configuration file.
#
#
# .. code-block:: json
#
#    {
#        "serializers": {
#            "ase.atoms.Atoms": "abc.ase.atoms.Atoms"
#        },
#        "deserializers": {
#            "aiida.orm.nodes.data.structure.StructureData": "aiida_pythonjob.data.deserializer.structure_data_to_pymatgen" # noqa
#        }
#    }
#
# The `orm.List`, `orm.Dict`and `orm.StructureData` data types already have built-in deserializers.
#

######################################################################
# What's Next
# -----------
# +-----------------------------------------+------------------------------------------------------+
# | `Tutorials <../tutorial/index.rst>`__   | Real-world examples in computational materials       |
# |                                         | science and more.                                    |
# |                                         |                                                      |
# +-----------------------------------------+------------------------------------------------------+
#
#
