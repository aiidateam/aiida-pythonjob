from aiida import load_profile
from aiida.engine import run
from aiida_pythonjob import PythonJob, prepare_pythonjob_inputs

load_profile()


def add(x, y):
    return x + y


inputs = prepare_pythonjob_inputs(
    add, function_inputs={"x": 1, "y": 2}, output_ports=[{"name": "add"}], computer="localhost"
)
run(PythonJob, inputs=inputs)
