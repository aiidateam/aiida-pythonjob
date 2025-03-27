import aiida
import pytest


def test_typing():
    """Test function with typing."""
    from typing import List

    from aiida_pythonjob.utils import get_required_imports
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
    from aiida_pythonjob.config import config
    from aiida_pythonjob.data.pickled_data import PickledData
    from aiida_pythonjob.data.serializer import serialize_to_aiida_nodes

    inputs = {"a": 1, "b": 2.0, "c": set()}
    with pytest.raises(
        ValueError,
        match="Cannot serialize type=set. No suitable method found",
    ):
        new_inputs = serialize_to_aiida_nodes(inputs)
    # Allow pickling
    config["allow_pickle"] = True
    new_inputs = serialize_to_aiida_nodes(inputs)
    assert isinstance(new_inputs["a"], aiida.orm.Int)
    assert isinstance(new_inputs["b"], aiida.orm.Float)
    assert isinstance(new_inputs["c"], PickledData)


def test_atoms_data():
    from aiida_pythonjob.data.atoms import AtomsData
    from ase.build import bulk

    atoms = bulk("Si")

    atoms_data = AtomsData(atoms)
    assert atoms_data.value == atoms


def test_only_data_with_value():
    from aiida_pythonjob.data import general_serializer

    # do not raise error because the built-in serializer can handle it
    general_serializer(aiida.orm.List([1]))
    # Test case: aiida.orm.ArrayData should raise a ValueError
    with pytest.raises(
        ValueError,
        match="AiiDA data: aiida.orm.nodes.data.array.xy.XyData, does not have a `value` attribute or deserializer.",
    ):
        general_serializer(aiida.orm.XyData())


def test_deserializer():
    import numpy as np
    from aiida_pythonjob.data.deserializer import deserialize_to_raw_python_data

    data = aiida.orm.ArrayData()
    data.set_array("data", np.array([1, 2, 3]))
    data = deserialize_to_raw_python_data(
        data,
        deserializers={
            "aiida.orm.nodes.data.array.array.ArrayData": "aiida_pythonjob.data.deserializer.generate_aiida_node_deserializer"  # noqa
        },
    )
    assert data == {"array|data": [3]}
