from aiida_pythonjob.utils import build_function_data


def test_build_function_data():
    from math import sqrt

    function_data = build_function_data(sqrt)
    assert function_data == {
        "name": "sqrt",
        "mode": "use_module_path",
        "source_code": "from math import sqrt",
    }
    #
    try:
        function_data = build_function_data(1)
    except Exception as e:
        assert str(e) == "Provided object is not a callable function or class."
