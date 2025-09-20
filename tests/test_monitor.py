import datetime

from aiida.engine import run_get_node

from aiida_pythonjob import MonitorPyFunction, prepare_monitor_function_inputs


def monitor_time(time: datetime.datetime):
    return datetime.datetime.now() > time


def test_async_function_runs_and_returns_result():
    inputs = prepare_monitor_function_inputs(
        monitor_time,
        function_inputs={"time": datetime.datetime.now() + datetime.timedelta(seconds=5)},
    )
    result, node = run_get_node(MonitorPyFunction, **inputs)
    assert node.is_finished_ok
    # The actual monitor function returns None
    assert result["result"].value is None


def test_async_function_raises_produces_exit_code():
    inputs = prepare_monitor_function_inputs(
        monitor_time,
        function_inputs={"time": datetime.datetime.now() + datetime.timedelta(seconds=20)},
        timeout=5.0,
    )
    _, node = run_get_node(MonitorPyFunction, **inputs)
    assert not node.is_finished_ok
    assert node.exit_status == MonitorPyFunction.exit_codes.ERROR_TIMEOUT.status
    assert "Monitor function execution timed out." in node.exit_message
