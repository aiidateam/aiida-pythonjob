from aiida.engine import run_get_node
from aiida_pythonjob import pyfunction


def test_function_default_outputs():
    """Test decorator."""

    @pyfunction()
    def add(x, y):
        return x + y

    result, node = run_get_node(add, x=1, y=2)

    assert result.value == 3
    assert node.process_label == "add"


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


def test_importable_function():
    """Test importable function."""
    from ase.build import bulk

    result, _ = run_get_node(pyfunction()(bulk), name="Si")
    assert result.value.get_chemical_formula() == "Si2"


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
            },
            {
                "name": "add_multiply.add",
                "identifier": "namespace",
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
        function_outputs=[
            {
                "name": "add_multiply",
                "identifier": "namespace",
            },
            {
                "name": "add_multiply.add",
                "identifier": "namespace",
            },
            {"name": "minus"},
        ],
    )

    assert result["add_multiply"]["add"]["order1"].value == 3
    assert result["add_multiply"]["add"]["order2"].value == 5
    assert result["add_multiply"]["multiply"].value == 2
