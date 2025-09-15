import numpy as np
import pytest
from aiida import orm

from aiida_pythonjob.data.serializer import all_serializers


class CustomData:
    """Data class that can be serialized to JSON."""

    def __init__(self, name, age, array):
        self.name = name
        self.age = age
        self.array = array

    def as_dict(self):
        return {"name": self.name, "age": self.age, "array": self.array}

    @classmethod
    def from_dict(cls, data):
        return cls(data["name"], data["age"], data["array"])


class NonJsonableData:
    """Data class that cannot be serialized to JSON."""

    def __init__(self, name, age):
        self.name = name
        self.age = age


@pytest.mark.parametrize(
    "data, data_type",
    (
        (1, orm.Int),
        (1.0, orm.Float),
        ("string", orm.Str),
        (np.float64(1.0), orm.Float),
        (np.int64(1), orm.Int),
        (np.bool_(True), orm.Bool),
        (True, orm.Bool),
        ([1, 2, 3], orm.List),
        ({"a": 1, "b": 2}, orm.Dict),
    ),
)
def test_serialize_aiida(data, data_type):
    from aiida_pythonjob.data.serializer import general_serializer

    serialized_data = general_serializer(data, serializers=all_serializers)
    assert isinstance(serialized_data, data_type)


def test_serialize_json():
    from aiida_pythonjob.data.jsonable_data import JsonableData
    from aiida_pythonjob.data.serializer import general_serializer

    data = CustomData("a", 1, np.zeros((3, 3)))

    serialized_data = general_serializer(data, serializers=all_serializers)
    assert isinstance(serialized_data, JsonableData)
