"""
PyFunction
===============

"""

######################################################################
# Default outputs
# --------------
#
# The default output of the function is `result`. The `pyfunction` task
# will store the result as one node in the database with the key `result`.
#
from aiida import load_profile
from aiida.engine import run_get_node
from aiida_pythonjob import pyfunction

load_profile()


@pyfunction()
def add(x, y):
    return x + y


result, node = run_get_node(add, x=1, y=2)
print("result: ", result)

######################################################################
# Custom outputs
# --------------
# If the function return a dictionary with fixed number of keys, and you
# want to store the values as separate outputs, you can specify the `outputs` parameter.
# For a dynamic number of outputs, you can use the namespace output, which is explained later.
#


@pyfunction(outputs=[{"name": "sum"}, {"name": "diff"}])
def add(x, y):
    return {"sum": x + y, "diff": x - y}


result, node = run_get_node(add, x=1, y=2)

print("result: ")
print("sum: ", result["sum"])
print("diff: ", result["diff"])

######################################################################
# Namespace Output
# --------------
#
# The `pyfunction` allows users to define namespace outputs. A namespace output
# is a dictionary with keys and values returned by a function. Each value in
# this dictionary will be serialized to AiiDA data, and the key-value pair
# will be stored in the database.
#
# Why Use Namespace Outputs?
#
# - **Dynamic and Flexible**: The keys and values in the namespace output are not fixed and can change based on the task's execution. # noqa
# - **Querying**: The data in the namespace output is stored as an AiiDA data node, allowing for easy querying and retrieval. # noqa
# - **Data Provenance**: When the data is used as input for subsequent tasks, the origin of data is tracked.
#
# For example: Consider a molecule adsorption calculation where the namespace
# output stores the surface slabs of the molecule adsorbed on different surface
# sites. The number of surface slabs can vary depending on the surface. These
# output surface slabs can be utilized as input to the next task to calculate the energy.

from ase import Atoms  # noqa: E402
from ase.build import bulk  # noqa: E402


@pyfunction(outputs=[{"name": "scaled_structures", "identifier": "namespace"}])
def generate_structures(structure: Atoms, factor_lst: list) -> dict:
    """Scale the structure by the given factor_lst."""
    scaled_structures = {}
    for i in range(len(factor_lst)):
        atoms = structure.copy()
        atoms.set_cell(atoms.cell * factor_lst[i], scale_atoms=True)
        scaled_structures[f"s_{i}"] = atoms
    return {"scaled_structures": scaled_structures}


result, node = run_get_node(generate_structures, structure=bulk("Al"), factor_lst=[0.95, 1.0, 1.05])
print("scaled_structures: ")
for key, value in result["scaled_structures"].items():
    print(key, value)


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


@pyfunction()
def add(x, y):
    sum = x + y
    if sum < 0:
        exit_code = {"status": 410, "message": "Some elements are negative"}
        return {"sum": sum, "exit_code": exit_code}
    return {"sum": sum}


result, node = run_get_node(add, x=1, y=-2)
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


@pyfunction()
def make_supercell(structure, n=2):
    return structure * [n, n, n]


structure = orm.StructureData(cell=[[1, 0, 0], [0, 1, 0], [0, 0, 1]])
structure.append_atom(position=(0.0, 0.0, 0.0), symbols="Li")

result, node = run_get_node(
    make_supercell,
    structure=structure,
    deserializers={
        "aiida.orm.nodes.data.structure.StructureData": "aiida_pythonjob.data.deserializer.structure_data_to_atoms"
    },
)
print("result: ", result)

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
