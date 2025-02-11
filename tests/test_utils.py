import pytest
from aiida_pythonjob.utils import build_function_data


def test_build_function_data():
    """Test the build_function_data function behavior."""

    with pytest.raises(TypeError, match="Provided object is not a callable function or class."):
        build_function_data(1)

    function_data = build_function_data(build_function_data)
    assert function_data["name"] == "build_function_data"
    assert "source_code" in function_data
    assert "pickled_function" in function_data
    node = function_data["pickled_function"]
    with node.base.repository.open(node.FILENAME, mode="rb") as f:
        text = f.read()
        assert b"cloudpickle" not in text

    function_data = build_function_data(build_function_data, register_pickle_by_value=True)
    assert function_data["name"] == "build_function_data"
    assert "source_code" in function_data
    assert "pickled_function" in function_data
    node = function_data["pickled_function"]
    with node.base.repository.open(node.FILENAME, mode="rb") as f:
        text = f.read()
        assert b"cloudpickle" in text

    def local_function(x, y):
        return x + y

    function_data = build_function_data(local_function)
    assert function_data["name"] == "local_function"
    assert "source_code" in function_data
    assert function_data["mode"] == "use_pickled_function"

    def outer_function():
        def nested_function(x, y):
            return x + y

        return nested_function

    nested_func = outer_function()
    function_data = build_function_data(nested_func)
    assert function_data["name"] == "nested_function"
    assert function_data["mode"] == "use_pickled_function"
