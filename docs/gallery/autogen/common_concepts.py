"""
=======================================================
Common concepts
=======================================================

.. _common_topics:

The ``PythonJob`` and ``pyfunction`` tools share a common API for several advanced features,
including input/output specification, data serialization, and custom exit
codes. This guide covers these shared concepts.

The examples below are illustratived using `PyFunction`, but the same principles
apply to `PythonJob` as well.
"""

# %%#
# Annotating inputs and outputs
# --------------------------------
# Both ``PythonJob`` and ``pyfunction`` use a specification system to define how
# function inputs and outputs are structured and stored as AiiDA nodes. This is
# essential for ensuring that returned data is properly serialized and can be
# queried or reused in subsequent workflows.
#
#
# Static namespaces
# -----------------
# For functions that return a dictionary with a fixed set of keys, you can define a static namespace.
# Each key-value pair in the returned dictionary will be stored as a separate output node.
#
#

from typing import Annotated, Any

from aiida import load_profile
from aiida.engine import run_get_node

from aiida_pythonjob import PyFunction, prepare_pyfunction_inputs, spec

load_profile()


def add_multiply(x, y):
    return {"sum": x + y, "product": x * y}


# Usage with PythonJob
inputs = prepare_pyfunction_inputs(
    add_multiply,
    function_inputs={"x": 1, "y": 2},
    outputs_spec=spec.namespace(sum=Any, product=Any),
)

result, node = run_get_node(PyFunction, inputs=inputs)
print("sum:", result["sum"])
print("product:", result["product"])


# %%
# One can also annotate the return type of the function using Python's type hints.
# Then the specification can be inferred automatically.
# For example, the function above can be defined as:


def add_multiply(x, y) -> Annotated[dict, spec.namespace(sum=Any, product=Any)]:
    return {"sum": x + y, "product": x * y}


inputs = prepare_pyfunction_inputs(
    add_multiply,
    function_inputs={"x": 1, "y": 2},
)

result, node = run_get_node(PyFunction, inputs=inputs)
print("sum:", result["sum"])
print("product:", result["product"])

# %%
# .. note::
#
#    For pyfunction, one can use pass the specification directly to the decorator:
#
#    .. code-block:: python
#
#        @pyfunction(outputs=spec.namespace(sum=Any, product=Any))
#        def add_multiply(x, y):
#            return {"sum": x + y, "product": x * y}
#    And then run it directly:
#
#    .. code-block:: python
#
#        from aiida.engine import run_get_node
#        result, node = run_get_node(add_multiply, x=1, y=2)
#
# One can annotate the inputs using Python's type hints as well:
#


def add_multiply(
    data: Annotated[dict, spec.namespace(x=int, y=int)],
) -> Annotated[dict, spec.namespace(sum=Any, product=Any)]:
    x = data["x"]
    y = data["y"]
    return {"sum": x + y, "product": x * y}


# Usage with PythonJob
inputs = prepare_pyfunction_inputs(
    add_multiply,
    function_inputs={"data": {"x": 1, "y": 2}},
    outputs_spec=spec.namespace(sum=Any, product=Any),
)

result, node = run_get_node(PyFunction, inputs=inputs)
print(node.inputs.function_inputs.data.x)
print(node.inputs.function_inputs.data.y)
print(node.outputs.sum)
print(node.outputs.product)

# %%
# Dynamic namespaces
# ~~~~~~~~~~~~~~~~~~~~~~
# When the number and names of outputs are not known until runtime, you can use a dynamic namespace.
# This is ideal for functions that generate a variable number of results.
#
#


def generate_square_numbers(n):
    """Generate a dict of square numbers."""
    return {f"square_{i}": i**2 for i in range(n)}


# Usage with PythonJob:
inputs = prepare_pyfunction_inputs(
    generate_square_numbers,
    function_inputs={"n": 5},
    outputs_spec=spec.dynamic(Any),
)

result, node = run_get_node(PyFunction, inputs=inputs)
print("result: ")
for key, value in result.items():
    print(f"{key}: {value}")

# %%
# Nested namespaces
# ~~~~~~~~~~~~~~~~~~~~~~
# Namespaces can be nested to represent complex, hierarchical data structures,
# allowing you to unpack nested dictionaries into a corresponding nested output structure.
#


def nested_dict_task(x, y):
    """Returns a nested dictionary with a corresponding nested namespace."""
    return {"sum": x + y, "nested": {"diff": x - y, "product": x * y}}


inputs = prepare_pyfunction_inputs(
    nested_dict_task,
    function_inputs={"x": 5, "y": 3},
    outputs_spec=spec.namespace(sum=int, nested=spec.namespace(diff=int, product=int)),
)

result, node = run_get_node(PyFunction, inputs=inputs)
print("result: ")
print("sum:", result["sum"])
print("nested diff:", result["nested"]["diff"])
print("nested product:", result["nested"]["product"])

# %%#
# Custom exit codes
# --------------------------------
# You can signal the status of a task by returning a special ``exit_code``
# dictionary from your function. If the status is non-zero, the process will be
# marked as failed with the corresponding exit status and message.
#
# .. code-block:: python
#
#     def check_sum(x, y):
#         total = x + y
#         if total < 0:
#             exit_code = {"status": 410, "message": "Sum is negative."}
#             return {"sum": total, "exit_code": exit_code}
#         return {"sum": total}
#
#     # This will result in a failed process with exit status 410
#     result, node = run_get_node(check_sum, x=1, y=-21)
#     print("exit_status:", node.exit_status) -> 410

# %%#
# Data serialization and deserialization
# ------------------------------------------------
#
# The system provides a flexible mechanism for serializing Python objects into
# AiiDA data nodes and deserializing AiiDA nodes back into Python objects.
#
# Automatic serialization
# ~~~~~~~~~~~~~~~~~~~~~~~
# When you provide standard Python objects as inputs, they are automatically
# serialized:
#
# 1. The system first searches for an AiiDA data entry point matching the
#    object's type (e.g., ``ase.atoms.Atoms``).
# 2. If no specific serializer is found, it attempts to store the data using
#    ``JsonableData``.
# 3. If the data is not JSON-serializable, it will be pickled using
#    ``PickledData`` **only if** ``use_pickle=True`` is set.
#
# Registering a custom serializer
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To support a custom data type (e.g., ``ase.atoms.Atoms``), you must register
# its corresponding AiiDA data class as an ``aiida.data`` entry point in your
# package's configuration (e.g., ``pyproject.toml``).
#
# .. code-block:: toml
#
#    [project.entry-points."aiida.data"]
#    myplugin.ase.atoms.Atoms = "myplugin.data.atoms:MyAtomsData"
#
# This entry point tells the system to use your ``MyAtomsData`` class to
# serialize any ``ase.atoms.Atoms`` object.
#
# Custom deserializers
# ~~~~~~~~~~~~~~~~~~~~
# When passing an AiiDA data node as an input, it should have a ``.value`` attribute
# that returns a simple Python type. For example, an ``Int`` node has a ``.value``
# attribute that returns a standard Python integer. Similarly, from aiida-core v2.7.0
# onwards, the ``Dict`` and ``List`` nodes also have ``.value`` attributes that return
# standard Python dictionaries and lists. The ``PickledData`` and ``JsonableData`` nodes also have a ``.value``
# attribute that returns the unpickled Python object.
#
# However, if the input AiiDA node does not have a ``.value`` attribute, you need to
# register a deserializer for it. This ensures the node can be converted into a
# Python object that the function can process.
#
# For example, you can pass the `deserializer` argument to
# ``prepare_pyfunction_inputs`` or ``prepare_pypythonjob_inputs``:
#
# .. code-block:: python
#
#     from aiida import orm
#
#     # This requires a deserializer to convert the StructureData node
#     # into an ASE Atoms object.
#     inputs = prepare_pyfunction_inputs(
#         make_supercell,
#         deserializers={
#             "aiida.orm.nodes.data.structure.StructureData": "aiida_pythonjob.data.deserializer.structure_data_to_atoms"  # noqa: E501
#         },
#         function_inputs={"structure": orm.load_node(PK)},
#     )
#
# Configuration via ``pythonjob.json``
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# You can manage serializers and deserializers globally by creating a
# ``pythonjob.json`` file in your AiiDA configuration directory (e.g., ``~/.aiida``).
#
# .. code-block:: json
#
#    {
#        "serializers": {
#            "ase.atoms.Atoms": "myplugin.ase.atoms.Atoms"
#        },
#        "deserializers": {
#            "aiida.orm.nodes.data.structure.StructureData": "aiida_pythonjob.data.deserializer.structure_data_to_pymatgen"  # noqa: E501
#        },
#        "use_pickle": false
#    }
#
# Using pickle for non-JSON Data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For one-off cases involving non-JSON-serializable data, you can enable
# pickle serialization by passing ``use_pickle=True``.
#
# .. code-block:: python
#
#     inputs = prepare_pyfunction_inputs(
#         my_function,
#         function_inputs={"data": MyCustomObject(1)},
#         use_pickle=True
#     )
#
