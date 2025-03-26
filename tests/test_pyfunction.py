from aiida import orm
from aiida.engine import run_get_node
from aiida_pythonjob import pyfunction


def test_function_default_outputs(fixture_localhost):
    """Test decorator."""

    @pyfunction()
    def add(x, y):
        return x + y

    result, node = run_get_node(add, x=1, y=2)

    assert result.value == 3
    assert node.process_label == "add"


def test_output_tuple():
    @pyfunction(
        outputs=[
            {"name": "sum"},
            {"name": "diff"},
        ]
    )
    def add(x, y):
        return x + y, x - y

    result, node = run_get_node(add, x=1, y=2)

    assert result["sum"].value == 3
    assert result["diff"].value == -1


def test_function_custom_outputs():
    """Test decorator."""

    @pyfunction(
        outputs=[
            {"name": "sum"},
            {"name": "diff"},
        ]
    )
    def add(x, y):
        return {"sum": x + y, "diff": x - y}

    result, node = run_get_node(add, x=1, y=2)

    assert result["sum"].value == 3
    assert result["diff"].value == -1
    assert node.process_label == "add"


def test_function_custom_inputs_outputs():
    @pyfunction(
        inputs=[{"name": "volumes", "identifier": "namespace"}, {"name": "energies", "identifier": "namespace"}],
        outputs=[{"name": "volumes", "identifier": "namespace"}, {"name": "energies", "identifier": "namespace"}],
    )
    def plot_eos(volumes, energies):
        return {"volumes": volumes, "energies": energies}

    _, node = run_get_node(plot_eos, volumes={"s_1": 1, "s_2": 2, "s_3": 3}, energies={"s_1": 1, "s_2": 2, "s_3": 3})
    assert node.inputs.function_inputs.volumes.s_1.value == 1
    assert node.outputs.volumes.s_1.value == 1


def test_importable_function():
    """Test importable function."""
    from ase.build import bulk

    result, _ = run_get_node(pyfunction()(bulk), name="Si", cubic=True)
    assert result.value.get_chemical_formula() == "Si8"


def test_kwargs_inputs():
    """Test function with kwargs."""

    @pyfunction(outputs=[{"name": "sum"}])
    def add(x, y=1, **kwargs):
        x += y
        for value in kwargs.values():
            x += value
        return x

    result, _ = run_get_node(add, x=1, y=2, a=3, b=4)
    assert result["sum"].value == 10


def test_namespace_output():
    """Test function with namespace output and input."""

    @pyfunction(
        outputs=[
            {
                "name": "add_multiply",
                "identifier": "namespace",
                "ports": [{"name": "add", "identifier": "namespace"}, "multiply"],
            },
            {"name": "minus"},
        ]
    )
    def myfunc(x, y):
        add = {"order1": x + y, "order2": x * x + y * y}
        return {
            "add_multiply": {"add": add, "multiply": x * y},
            "minus": x - y,
        }

    result, node = run_get_node(myfunc, x=1, y=2)
    print("result: ", result)

    assert result["add_multiply"]["add"]["order1"].value == 3
    assert result["add_multiply"]["add"]["order2"].value == 5
    assert result["add_multiply"]["multiply"].value == 2


def test_override_outputs():
    """Test function with namespace output and input."""

    @pyfunction()
    def myfunc(x, y):
        add = {"order1": x + y, "order2": x * x + y * y}
        return {
            "add_multiply": {"add": add, "multiply": x * y},
            "minus": x - y,
        }

    result, node = run_get_node(
        myfunc,
        x=1,
        y=2,
        output_ports=[
            {
                "name": "add_multiply",
                "identifier": "namespace",
                "ports": [{"name": "add", "identifier": "namespace"}],
            },
            {"name": "minus"},
        ],
    )

    assert result["add_multiply"]["add"]["order1"].value == 3
    assert result["add_multiply"]["add"]["order2"].value == 5
    assert result["add_multiply"]["multiply"].value == 2


def test_function_execution_failed():
    @pyfunction()
    def add(x):
        import math

        return math.sqrt(x)

    _, node = run_get_node(add, x=-2)
    assert node.exit_status == 325


def test_exit_code():
    """Test function with exit code."""
    from numpy import array

    @pyfunction()
    def add(x: array, y: array) -> array:
        sum = x + y
        if (sum < 0).any():
            exit_code = {"status": 410, "message": "Some elements are negative"}
            return {"sum": sum, "exit_code": exit_code}
        return {"sum": sum}

    result, node = run_get_node(add, x=array([1, 1]), y=array([1, -2]))
    assert node.exit_status == 410
    assert node.exit_message == "Some elements are negative"


def test_aiida_node_as_inputs_outputs():
    """Test function with AiiDA nodes as inputs and outputs."""

    @pyfunction()
    def add(x, y):
        return {"sum": orm.Int(x + y), "diff": orm.Int(x - y)}

    result, node = run_get_node(add, x=orm.Int(1), y=orm.Int(2))
    assert set(result.keys()) == {"sum", "diff"}
    assert result["sum"].value == 3


def test_missing_output():
    @pyfunction(
        outputs=[
            {"name": "sum"},
            {"name": "diff"},
        ]
    )
    def add(x, y):
        return {"sum": x + y}

    result, node = run_get_node(add, x=1, y=2)

    assert node.exit_status == 11


def test_nested_inputs_outputs():
    """Test function with nested inputs and outputs."""

    @pyfunction(
        inputs=[
            {
                "name": "input1",
                "identifier": "namespace",
                "ports": [
                    {"name": "x1"},
                    {"name": "y1"},
                ],
            },
            {
                "name": "input1",
                "identifier": "namespace",
                "ports": [
                    {"name": "x2"},
                    {"name": "y2"},
                ],
            },
        ],
        outputs=[
            {
                "name": "result1",
                "identifier": "namespace",
                "ports": [
                    {"name": "sum1"},
                    {"name": "diff1"},
                ],
            },
            {
                "name": "result2",
                "identifier": "namespace",
                "ports": [
                    {"name": "sum2"},
                    {"name": "diff2"},
                ],
            },
        ],
    )
    def add(input1, input2):
        return {
            "result1": {"sum1": input1["x"] + input1["y"], "diff1": input1["x"] - input1["y"]},
            "result2": {"sum2": input2["x"] + input2["y"], "diff2": input2["x"] - input2["y"]},
        }

    result, node = run_get_node(add, input1={"x": 1, "y": 2}, input2={"x": 1, "y": 3})

    assert node.outputs.result1.sum1.value == 3
    assert node.outputs.result1.diff1.value == -1
    assert node.outputs.result2.sum2.value == 4
    assert node.outputs.result2.diff2.value == -2
