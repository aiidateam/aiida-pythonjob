import os
import pathlib
import tempfile

import pytest
from aiida import orm
from aiida.engine import run_get_node
from aiida_pythonjob import PythonJob, prepare_pythonjob_inputs


def test_validate_inputs(fixture_localhost):
    def add(x, y):
        return x + y

    with pytest.raises(ValueError, match="Either function or function_data must be provided"):
        prepare_pythonjob_inputs(
            function_inputs={"x": 1, "y": 2},
        )
    with pytest.raises(ValueError, match="Only one of function or function_data should be provided"):
        prepare_pythonjob_inputs(
            function=add,
            function_data={"module_path": "math", "name": "sqrt", "is_pickle": False},
        )


def test_function_default_outputs(fixture_localhost):
    """Test decorator."""

    def add(x, y):
        return x + y

    inputs = prepare_pythonjob_inputs(
        add,
        function_inputs={"x": 1, "y": 2},
        process_label="add",
    )
    result, node = run_get_node(PythonJob, **inputs)

    assert result["result"].value == 3
    assert node.process_label == "add"


def test_function_custom_outputs(fixture_localhost):
    """Test decorator."""

    def add(x, y):
        return {"sum": x + y, "diff": x - y}

    inputs = prepare_pythonjob_inputs(
        add,
        function_inputs={"x": 1, "y": 2},
        output_ports=[
            {"name": "sum"},
            {"name": "diff"},
        ],
    )
    result, node = run_get_node(PythonJob, **inputs)

    assert result["sum"].value == 3
    assert result["diff"].value == -1
    assert node.process_label == "PythonJob<add>"


@pytest.mark.skip("Can not inspect the built-in function.")
def test_importable_function(fixture_localhost):
    """Test importable function."""
    from operator import add

    inputs = prepare_pythonjob_inputs(
        add,
        function_inputs={"x": 1, "y": 2},
        output_ports=[
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
        output_ports=[
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
        output_ports=[
            {
                "name": "add_multiply",
                "identifier": "namespace",
                "ports": [{"name": "add", "identifier": "namespace"}, "multiply"],
            },
            {"name": "minus"},
        ],
    )
    result, node = run_get_node(PythonJob, **inputs)
    print("result: ", result)

    assert result["add_multiply"]["add"]["order1"].value == 3
    assert result["add_multiply"]["add"]["order2"].value == 5
    assert result["add_multiply"]["multiply"].value == 2


def test_parent_folder_remote(fixture_localhost):
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
        output_ports=[{"name": "sum"}],
    )
    result1, node1 = run_get_node(PythonJob, inputs=inputs1)

    inputs2 = prepare_pythonjob_inputs(
        multiply,
        function_inputs={"x": 1, "y": 2},
        output_ports=[{"name": "product"}],
        parent_folder=result1["remote_folder"],
    )
    result2, node2 = run_get_node(PythonJob, inputs=inputs2)

    assert result2["product"].value == 5


def test_parent_folder_local(fixture_localhost):
    """Test function with parent folder."""

    with tempfile.TemporaryDirectory() as tmpdir:
        dirpath = pathlib.Path(tmpdir)
        with open((dirpath / "result.txt"), "w") as f:
            f.write("3")

        parent_folder = orm.FolderData(tree=dirpath.absolute())

        def multiply(x, y):
            with open("parent_folder/result.txt", "r") as f:
                z = int(f.read())
            return x * y + z

        inputs2 = prepare_pythonjob_inputs(
            multiply,
            function_inputs={"x": 1, "y": 2},
            output_ports=[{"name": "product"}],
            parent_folder=parent_folder,
        )
        result2, node2 = run_get_node(PythonJob, inputs=inputs2)

        assert result2["product"].value == 5


def test_upload_files(fixture_localhost):
    """Test function with upload files."""

    def add():
        with open("input.txt", "r") as f:
            a = int(f.read())
        with open("another_input.txt", "r") as f:
            b = int(f.read())
        with open("inputs_folder/another_input.txt", "r") as f:
            c = int(f.read())
        return a + b + c

    # create a temporary file "input.txt" in the current directory
    with tempfile.TemporaryDirectory() as tmpdir:
        dirpath = pathlib.Path(tmpdir)
        with open((dirpath / "input.txt"), "w") as f:
            f.write("2")
        with open((dirpath / "another_input.txt"), "w") as f:
            f.write("3")
        # create a temporary folder "inputs_folder"
        os.makedirs((dirpath / "inputs_folder"), exist_ok=True)
        with open((dirpath / "inputs_folder/another_input.txt"), "w") as f:
            f.write("4")
        # we need use full path to the file
        input_file = str(dirpath / "input.txt")
        input_folder = str(dirpath / "inputs_folder")
        single_file_data = orm.SinglefileData(file=(dirpath / "another_input.txt"))
        # ------------------------- Submit the calculation -------------------
        inputs = prepare_pythonjob_inputs(
            add,
            upload_files={
                "input.txt": input_file,
                "another_input.txt": single_file_data,
                "inputs_folder": input_folder,
            },
        )
        result, node = run_get_node(PythonJob, inputs=inputs)
        assert result["result"].value == 9


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


def test_copy_files(fixture_localhost):
    """Test function with copy files."""

    def add(x, y):
        z = x + y
        with open("result.txt", "w") as f:
            f.write(str(z))

    def multiply(x_folder_name, y):
        with open(f"{x_folder_name}/result.txt", "r") as f:
            x = int(f.read())
        return x * y

    inputs = prepare_pythonjob_inputs(add, function_inputs={"x": 1, "y": 2})
    result, node = run_get_node(PythonJob, inputs=inputs)
    inputs = prepare_pythonjob_inputs(
        multiply,
        function_inputs={"x_folder_name": "x_folder_name", "y": 2},
        copy_files={"x_folder_name": result["remote_folder"]},
    )
    result, node = run_get_node(PythonJob, inputs=inputs)
    assert result["result"].value == 6


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


def test_local_function(fixture_localhost):
    def multily(x, y):
        return x * y

    def add(x, y):
        return x + multily(x, y)

    inputs = prepare_pythonjob_inputs(
        add,
        function_inputs={"x": 2, "y": 3},
    )
    result, node = run_get_node(PythonJob, **inputs)
    assert result["result"].value == 8


@pytest.mark.usefixtures("started_daemon_client")
def test_submit(fixture_localhost):
    """Test decorator."""
    from aiida.engine import submit

    def add(x, y):
        return x + y

    inputs = prepare_pythonjob_inputs(
        add,
        function_inputs={"x": 1, "y": 2},
        process_label="add",
    )
    node = submit(PythonJob, **inputs, wait=True)

    assert node.outputs.result.value == 3
    assert node.process_label == "add"
