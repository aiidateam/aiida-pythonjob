import aiida
from aiida_pythonjob.data import general_serializer
from aiida_pythonjob.utils import get_required_imports


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

    modules = get_required_imports(generate_structures)
    assert modules == {
        "typing": {"List"},
        "builtins": {"list", "float"},
        "numpy": {"array"},
    }


def test_python_job():
    """Test a simple python node."""
    from aiida_pythonjob.data.pickled_data import PickledData
    from aiida_pythonjob.data.serializer import serialize_to_aiida_nodes

    inputs = {"a": 1, "b": 2.0, "c": set()}
    new_inputs = serialize_to_aiida_nodes(inputs)
    assert isinstance(new_inputs["a"], aiida.orm.Int)
    assert isinstance(new_inputs["b"], aiida.orm.Float)
    assert isinstance(new_inputs["c"], PickledData)


def test_dict_list():
    from aiida_pythonjob.data.data_with_value import Dict, List

    data = List([1, 2, 3])
    assert data.value == [1, 2, 3]
    data = Dict({"a": 1, "b": 2})
    assert data.value == {"a": 1, "b": 2}


def test_atoms_data():
    from aiida_pythonjob.data.atoms import AtomsData
    from ase.build import bulk

    atoms = bulk("Si")

    atoms_data = AtomsData(atoms)
    assert atoms_data.value == atoms


def test_only_data_with_value():
    try:
        general_serializer(aiida.orm.List([1]))
    except ValueError as e:
        assert str(e) == "Only AiiDA data Node with a value attribute is allowed."
