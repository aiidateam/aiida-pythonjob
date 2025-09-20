import asyncio
import time

import pytest
from aiida.engine import run_get_node, submit
from aiida.engine.processes import control

from aiida_pythonjob import PyFunction, prepare_pyfunction_inputs, pyfunction


@pyfunction()
async def add_async(x, y, t=0.01):
    await asyncio.sleep(t)
    return x + y


@pyfunction()
async def fail_async(x):
    await asyncio.sleep(0)
    raise ValueError(f"bad x: {x}")


def test_async_function_runs_and_returns_result():
    inputs = prepare_pyfunction_inputs(
        add_async,
        function_inputs={"x": 1, "y": 2},
    )
    result, node = run_get_node(PyFunction, **inputs)
    assert node.is_finished_ok
    assert "result" in result
    assert result["result"].value == 3


def test_async_function_raises_produces_exit_code():
    inputs = prepare_pyfunction_inputs(
        fail_async,
        function_inputs={"x": 99},
    )
    _, node = run_get_node(PyFunction, **inputs)
    assert not node.is_finished_ok
    assert node.exit_status == PyFunction.exit_codes.ERROR_FUNCTION_EXECUTION_FAILED.status
    assert "Function execution failed." in node.exit_message
    assert "bad x: 99" in node.exit_message


@pytest.mark.usefixtures("started_daemon_client")
def test_async_function_kill():
    inputs = prepare_pyfunction_inputs(
        add_async,
        function_inputs={"x": 1, "y": 2, "t": 30},
    )
    node = submit(PyFunction, **inputs)
    # wait process to start
    time.sleep(5)
    control.kill_processes(
        [node],
        all_entries=None,
        timeout=10,
    )
    # wait kill to take effect
    time.sleep(10)
    assert node.is_killed
    assert "Killed through" in node.process_status
