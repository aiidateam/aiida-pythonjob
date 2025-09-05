"""
PyFunction
===============

"""

# %%#
# Default and custom outputs
# --------------------------
#
# The ``@pyfunction`` decorator turns a Python function into an AiiDA ``PyFunction``.
# By default, the entire return value is stored in a single output node named ``result``.
#
from aiida import load_profile
from aiida.engine import run_get_node

from aiida_pythonjob import pyfunction, spec

load_profile()


@pyfunction()
def add_and_subtract(x, y):
    return {"sum": x + y, "diff": x - y}


result, node = run_get_node(add_and_subtract, x=1, y=2)
print("Default result: ", result)


# %%#
# You can specify the ``outputs`` parameter to unpack a returned dictionary into
# separate output nodes.
#
@pyfunction(outputs=spec.namespace(sum=any, diff=any))
def add_and_subtract(x, y):
    return {"sum": x + y, "diff": x - y}


result, node = run_get_node(add_and_subtract, x=1, y=2)
print("Unpacked results: ")
print("sum: ", result["sum"])
print("diff: ", result["diff"])

# %%#
# Advanced features
# -----------------
# ``pyfunction`` supports many advanced features for data handling and workflow control.
# These functionalities are shared with ``PythonJob``.
#
# .. seealso::
#
#    For a detailed guide on the following topics, please refer to the :doc:`Common Concepts <common_concepts>`:
#
#    - **Dynamic and Nested Namespaces**: For handling complex or variable outputs.
#    - **Custom Exit Codes**: For robust error handling and workflow control.
#    - **Data Serialization and Deserialization**: For working with custom data types like ``ase.Atoms`` or other AiiDA data nodes.  # noqa: E501
#
# Here is a more complex example demonstrating a dynamic namespace output with custom data types.
from ase import Atoms  # noqa: E402
from ase.build import bulk  # noqa: E402


@pyfunction(outputs=spec.dynamic(Atoms))
def generate_structures(element: str, factors: list) -> dict:
    """Scale a bulk structure by the given factors."""
    scaled_structures = {}
    initial_structure = bulk(element)
    for i, factor in enumerate(factors):
        atoms = initial_structure.copy()
        atoms.set_cell(atoms.cell * factor, scale_atoms=True)
        scaled_structures[f"s_{i}"] = atoms
    return scaled_structures


result, node = run_get_node(generate_structures, element="Al", factors=[0.95, 1.0, 1.05])

print("Generated scaled structures:")
for key, value in result.items():
    print(key, value)
