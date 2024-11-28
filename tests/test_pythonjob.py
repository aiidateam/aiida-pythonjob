import pytest
from aiida.engine import run_get_node
from aiida_pythonjob import PythonJob, prepare_pythonjob_inputs


def test_function_default_outputs(fixture_localhost):
    """Test decorator."""

    def add(x, y):
        return x + y

    inputs = prepare_pythonjob_inputs(
        add,
        function_inputs={"x": 1, "y": 2},
    )
    result, node = run_get_node(PythonJob, **inputs)
    print("result: ", result)

    assert result["result"].value == 3
    assert node.process_label == "PythonJob<add>"


def test_function_custom_outputs(fixture_localhost):
    """Test decorator."""

    def add(x, y):
        return {"sum": x + y, "diff": x - y}

    inputs = prepare_pythonjob_inputs(
        add,
        function_inputs={"x": 1, "y": 2},
        function_outputs=[
            {"name": "sum"},
            {"name": "diff"},
        ],
    )
    result, node = run_get_node(PythonJob, **inputs)

    assert result["sum"].value == 3
    assert result["diff"].value == -1


@pytest.mark.skip("Can not inspect the built-in function.")
def test_importable_function(fixture_localhost):
    """Test importable function."""
    from operator import add

    inputs = prepare_pythonjob_inputs(
        add,
        function_inputs={"x": 1, "y": 2},
        function_outputs=[
            {"name": "sum"},
        ],
    )
    result, node = run_get_node(PythonJob, **inputs)
    print("result: ", result)
    assert result["sum"].value == 3


def test_kwargs_inputs(fixture_localhost):
    """Test function with kwargs."""

    def add(x, y=1, **kwargs):
        x += y
        for value in kwargs.values():
            x += value
        return x

    inputs = prepare_pythonjob_inputs(
        add,
        function_inputs={"x": 1, "y": 2, "a": 3, "b": 4},
        function_outputs=[
            {"name": "sum"},
        ],
    )
    result, node = run_get_node(PythonJob, **inputs)
    assert result["sum"].value == 10


def test_namespace_output(fixture_localhost):
    """Test function with namespace output and input."""

    def myfunc(x, y):
        add = {"order1": x + y, "order2": x * x + y * y}
        return {
            "add_multiply": {"add": add, "multiply": x * y},
            "minus": x - y,
        }

    inputs = prepare_pythonjob_inputs(
        myfunc,
        function_inputs={"x": 1, "y": 2},
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
    result, node = run_get_node(PythonJob, **inputs)
    print("result: ", result)

    assert result["add_multiply"]["add"]["order1"].value == 3
    assert result["add_multiply"]["add"]["order2"].value == 5
    assert result["add_multiply"]["multiply"].value == 2


def test_parent_folder(fixture_localhost):
    """Test function with parent folder."""

    def add(x, y):
        z = x + y
        with open("result.txt", "w") as f:
            f.write(str(z))
        return x + y

    def multiply(x, y):
        with open("parent_folder/result.txt", "r") as f:
            z = int(f.read())
        return x * y + z

    inputs1 = prepare_pythonjob_inputs(
        add,
        function_inputs={"x": 1, "y": 2},
        function_outputs=[{"name": "sum"}],
    )
    result1, node1 = run_get_node(PythonJob, inputs=inputs1)

    inputs2 = prepare_pythonjob_inputs(
        multiply,
        function_inputs={"x": 1, "y": 2},
        function_outputs=[{"name": "product"}],
        parent_folder=result1["remote_folder"],
    )
    result2, node2 = run_get_node(PythonJob, inputs=inputs2)

    assert result2["product"].value == 5


def test_upload_files(fixture_localhost):
    """Test function with upload files."""

    # create a temporary file "input.txt" in the current directory
    with open("input.txt", "w") as f:
        f.write("2")

    # create a temporary folder "inputs_folder" in the current directory
    # and add a file "another_input.txt" in the folder
    import os

    os.makedirs("inputs_folder", exist_ok=True)
    with open("inputs_folder/another_input.txt", "w") as f:
        f.write("3")

    def add():
        with open("input.txt", "r") as f:
            a = int(f.read())
        with open("inputs_folder/another_input.txt", "r") as f:
            b = int(f.read())
        return a + b

    # ------------------------- Submit the calculation -------------------
    # we need use full path to the file
    input_file = os.path.abspath("input.txt")
    input_folder = os.path.abspath("inputs_folder")
    inputs = prepare_pythonjob_inputs(
        add,
        upload_files={
            "input.txt": input_file,
            "inputs_folder": input_folder,
        },
    )
    result, node = run_get_node(PythonJob, inputs=inputs)

    # wait=True)
    assert result["result"].value == 5


def test_retrieve_files(fixture_localhost):
    """Test retrieve files."""

    def add(x, y):
        z = x + y
        with open("result.txt", "w") as f:
            f.write(str(z))
        return x + y

    inputs = prepare_pythonjob_inputs(
        add,
        function_inputs={"x": 1, "y": 2},
        metadata={
            "options": {
                "additional_retrieve_list": ["result.txt"],
            }
        },
    )
    result, node = run_get_node(PythonJob, inputs=inputs)
    # ------------------------- Submit the calculation -------------------

    assert "result.txt" in result["retrieved"].list_object_names()


def test_exit_code(fixture_localhost):
    """Test function with exit code."""
    from numpy import array

    def add(x: array, y: array) -> array:
        sum = x + y
        if (sum < 0).any():
            exit_code = {"status": 410, "message": "Some elements are negative"}
            return {"sum": sum, "exit_code": exit_code}
        return {"sum": sum}

    inputs = prepare_pythonjob_inputs(
        add,
        function_inputs={"x": array([1, 1]), "y": array([1, -2])},
    )
    result, node = run_get_node(PythonJob, inputs=inputs)
    assert node.exit_status == 410
    assert node.exit_message == "Some elements are negative"
