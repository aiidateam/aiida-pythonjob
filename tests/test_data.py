from aiida_pythonjob import PickledFunction


def test_typing():
    """Test function with typing."""
    from typing import List

    from numpy import array

    def generate_structures(
        strain_lst: List[float],
        data: array,
        data1: array,
        strain_lst1: list,
    ) -> list[array]:
        pass

    modules = PickledFunction.get_required_imports(generate_structures)
    assert modules == {
        "typing": {"List"},
        "builtins": {"list", "float"},
        "numpy": {"array"},
    }
