import pytest
from aiida import orm
from aiida.engine import run_get_node, run_get_pk
from node_graph import socket_spec as spec

from aiida_pythonjob import PyFunction, prepare_pyfunction_inputs, pyfunction


@pyfunction()
def add(x, y):
    return x + y


def test_function_default_outputs(fixture_localhost):
    """Test decorator."""

    result, node = run_get_node(add, x=1, y=2)

    assert result.value == 3
    assert node.process_label == "add"


def test_run_get_pk(fixture_localhost):
    """Test decorator."""

    result, pk = run_get_pk(add, x=1, y=2)

    assert result.value == 3
    assert isinstance(pk, int)


def test_function_none_outputs(fixture_localhost):
    @pyfunction()
    def add(x, y):
        """Does not return anything."""

    result, node = run_get_node(add, x=1, y=2)
    assert result == {}
    assert node.is_finished_ok


def test_prepare_pyfunction_inputs():
    """Test prepare_pyfunction_inputs utility function."""
    inputs = prepare_pyfunction_inputs(
        add,
        function_inputs={"x": 1, "y": 2},
    )
    result, _ = run_get_node(PyFunction, **inputs)
    assert result["result"].value == 3


def test_output_tuple():
    @pyfunction(outputs=spec.namespace(sum=int, diff=int))
    def add(x, y):
        return x + y, x - y

    result, _ = run_get_node(add, x=1, y=2)

    assert result["sum"].value == 3
    assert result["diff"].value == -1


def test_function_custom_outputs():
    """Test decorator."""

    @pyfunction(outputs=spec.namespace(sum=int, diff=int))
    def add(x, y):
        return {"sum": x + y, "diff": x - y}

    result, node = run_get_node(add, x=1, y=2)

    assert result["sum"].value == 3
    assert result["diff"].value == -1
    assert node.process_label == "add"


def test_function_custom_inputs_outputs():
    @pyfunction(
        inputs=spec.namespace(volumes=spec.dynamic(any), energies=spec.dynamic(any)),
        outputs=spec.namespace(volumes=spec.dynamic(any), energies=spec.dynamic(any)),
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

    @pyfunction()
    def add(x, y=1, **kwargs):
        x += y
        for value in kwargs.values():
            x += value
        return x

    result, _ = run_get_node(add, x=1, y=2, a=3, b=4)
    assert result.value == 10


def test_namespace_output():
    """Test function with namespace output and input."""

    out = spec.namespace(
        add_multiply=spec.namespace(add=spec.dynamic(any), multiply=any),
        minus=any,
    )

    @pyfunction(outputs=out)
    def myfunc(x, y):
        """Function that returns a namespace output."""
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
        outputs_spec=spec.namespace(
            add_multiply=spec.namespace(add=spec.dynamic(any), multiply=any),
            minus=any,
        ),
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
    assert node.exit_status == 323


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
    print("result: ", result)
    assert set(result.keys()) == {"sum", "diff"}
    assert result["sum"].value == 3


def test_missing_output():
    @pyfunction(outputs=spec.namespace(sum=int, diff=int))
    def add(x, y):
        return {"sum": x + y}

    result, node = run_get_node(add, x=1, y=2)

    assert node.exit_status == 11


def test_nested_inputs_outputs():
    """Test function with nested inputs and outputs."""

    inp = spec.namespace(
        input1=spec.namespace(x=int, y=int),
        input2=spec.namespace(x=int, y=int),
    )
    out = spec.namespace(
        result1=spec.namespace(sum1=int, diff1=int),
        result2=spec.namespace(sum2=int, diff2=int),
    )

    @pyfunction(inputs=inp, outputs=out)
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


def test_top_level_outputs_dynamic():
    """Test function with dynamic top-level outputs."""

    @pyfunction(outputs=spec.dynamic(any))
    def test_dynamic(n: int):
        return {f"data_{i}": i for i in range(n)}

    result, node = run_get_node(
        test_dynamic,
        n=3,
    )

    # outputs should be serialized as dynamic rows
    assert node.outputs.data_0.value == 0
    assert node.outputs.data_1.value == 1


def test_dynamic_rows():
    """Test function with dynamic rows."""

    row = spec.namespace(sum=any, product=any)

    @pyfunction(outputs=spec.dynamic(row, sum=int))
    def test_dynamic_rows(data: spec.dynamic(row, sum=int)):
        return data

    result, node = run_get_node(
        test_dynamic_rows,
        data={
            "data_0": {"sum": 0, "product": 0},
            "data_1": {"sum": 1, "product": 2},
            "data_2": {"sum": 2, "product": 4},
            "sum": 1,
        },
    )

    # inputs should be serialized as dynamic rows
    assert node.inputs.function_inputs.data.sum.value == 1
    assert node.inputs.function_inputs.data.data_0.sum.value == 0
    # outputs should be serialized as dynamic rows
    assert node.outputs.sum.value == 1
    assert node.outputs.data_0.sum.value == 0


def test_only_data_with_value():
    from aiida import orm

    @pyfunction()
    def add(x, y):
        return x + y

    with pytest.raises(ValueError, match="Cannot deserialize AiiDA data of type"):
        run_get_node(add, x=1, y=orm.XyData())
