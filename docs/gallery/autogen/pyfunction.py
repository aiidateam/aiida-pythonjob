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


# %%
# Async functions
# ---------------
# ``pyfunction`` also supports Python's ``async`` functions. This is a powerful feature for
# tasks that are I/O-bound (e.g., waiting for network requests, file operations) or for
# running multiple tasks concurrently without blocking the AiiDA daemon.
#
# When you ``submit`` an async function, the call returns immediately with a process node,
# allowing your script to continue running while the function executes in the background.
#

from aiida.engine import submit
import datetime
from aiida_pythonjob import prepare_pyfunction_inputs


@pyfunction()
async def add_async(x, y, time: float):
    """A simple function that adds two numbers."""
    import asyncio

    # Simulate asynchronous I/O or computation
    await asyncio.sleep(time)
    return x + y


inputs = prepare_pyfunction_inputs(
    add_async,
    function_inputs={"x": 2, "y": 3, "time": 2.0},
)

node = submit(add_async, **inputs)

# %%
# Here is an example to monitor external events or conditions without blocking.
# Here is an example that waits until a specified time.
#


@pyfunction()
async def monitor_time(time: datetime.datetime):
    import asyncio

    # monitor until the specified time
    while datetime.datetime.now() < time:
        print("Waiting...")
        await asyncio.sleep(0.5)


inputs = prepare_pyfunction_inputs(
    monitor_time,
    function_inputs={"time": datetime.datetime.now() + datetime.timedelta(seconds=5)},
)

node = submit(monitor_time, **inputs)

# %%#
# Killing an async process
# ------------------------
# Since async functions run as regular AiiDA processes, they can be controlled and killed
# programmatically. This is useful for managing long-running or stuck tasks.
# You can kill a running async function using the AiiDA command line interface.
#
# .. code-block:: bash
#
#    $ verdi process kill <pk>
#
