import asyncio

from aiida.engine import run_get_node

from aiida_pythonjob import PyFunction, prepare_pyfunction_inputs, pyfunction


@pyfunction()
async def add_async(x, y):
    await asyncio.sleep(0.01)
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
