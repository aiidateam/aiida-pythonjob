import math

import numpy as np
import pytest

from aiida_pythonjob.data.jsonable_data import JsonableData


class WithToDictFromDict:
    """Uses to_dict / from_dict (classmethod) pair; contains diverse data types."""

    def __init__(self, a: int, arr, scal, floats):
        self.a = a
        self.arr = np.array(arr)
        self.scal = np.float32(scal)
        self.floats = floats  # dict with NaN/Inf

    def to_dict(self):
        return {
            "a": self.a,
            "arr": self.arr,  # ndarray
            "scal": self.scal,  # np.generic scalar
            "floats": self.floats,  # includes NaN/Inf
            "nested": {"tup": (1, 2)},  # tuple -> will serialize as list
        }

    @classmethod
    def from_dict(cls, d):
        # Accept list for "arr" and "nested.tup" after JSON round-trip
        arr = np.array(d["arr"])
        scal = np.float32(d["scal"])
        # "tup" may come back as list; normalize to tuple
        if isinstance(d["nested"]["tup"], list):
            d["nested"]["tup"] = tuple(d["nested"]["tup"])
        obj = cls(d["a"], arr, scal, d["floats"])
        return obj


class WithAsDictFromdict:
    """Uses as_dict / fromdict (alt names)."""

    def __init__(self, x: str):
        self.x = x

    def as_dict(self):
        return {"x": self.x}

    @classmethod
    def fromdict(cls, d):
        return cls(d["x"])


class MissingFromDict:
    """Has to_dict but *no* from_dict/fromdict -> should fail validation."""

    def to_dict(self):
        return {"y": 42}


def test_jsonabledata_stores_attributes_and_roundtrips_numpy_and_float_constants():
    obj = WithToDictFromDict(
        a=7,
        arr=[[1, 2], [3, 4]],
        scal=np.float32(5.5),
        floats={"pos": float("inf"), "neg": -float("inf"), "nan": float("nan")},
    )

    node = JsonableData(obj)

    attrs = node.base.attributes.all
    # Class / module metadata present
    assert "@class" in attrs and "@module" in attrs
    assert attrs["@class"] == "WithToDictFromDict"

    # NumPy array serialized to JSON-friendly form
    assert attrs["arr"] == [[1, 2], [3, 4]]
    # NumPy scalar converted to plain Python number
    assert isinstance(attrs["scal"], (int, float))

    # Special float tokens survive round-trip as strings before reconstruction
    # After json.dumps/loads(parse_constant=...), they should be strings
    assert attrs["floats"]["pos"] in ("Infinity", float("inf"))
    assert attrs["floats"]["neg"] in ("-Infinity", -float("inf"))
    # NaN becomes "NaN" after loads(parse_constant=lambda x: x)
    assert attrs["floats"]["nan"] in ("NaN", float("nan"))

    # Simulate reload (force reconstruction path)
    delattr(node, "_obj")
    rebuilt = node.obj
    assert isinstance(rebuilt, WithToDictFromDict)
    # Array restored as ndarray with equal contents
    np.testing.assert_array_equal(rebuilt.arr, np.array([[1, 2], [3, 4]]))
    # np.float32 preserved semantically
    assert math.isclose(float(rebuilt.scal), 5.5, rel_tol=1e-7)
    # Special floats correctly restored
    assert math.isinf(rebuilt.floats["pos"]) and rebuilt.floats["pos"] > 0
    assert math.isinf(rebuilt.floats["neg"]) and rebuilt.floats["neg"] < 0
    assert math.isnan(rebuilt.floats["nan"])


def test_supports_alternative_method_names_asdict_fromdict():
    obj = WithAsDictFromdict("hello")
    node = JsonableData(obj)
    delattr(node, "_obj")
    rebuilt = node.obj
    assert isinstance(rebuilt, WithAsDictFromdict)
    assert rebuilt.x == "hello"


def test_validation_requires_from_dict_like_method():
    obj = MissingFromDict()
    with pytest.raises(ValueError, match="The object must have at least one of the following methods"):
        JsonableData(obj)


def test_import_error_when_module_cannot_be_imported(monkeypatch):
    obj = WithAsDictFromdict("x")
    node = JsonableData(obj)
    attrs = node.base.attributes
    attrs.set("@module", "nonexistent.module.path")
    if hasattr(node, "_obj"):
        delattr(node, "_obj")
    with pytest.raises(ImportError):
        _ = node.obj


def test_attribute_error_when_class_missing(monkeypatch):
    obj = WithAsDictFromdict("x")
    node = JsonableData(obj)
    attrs = node.base.attributes
    attrs.set("@class", "NotAClassHere")
    if hasattr(node, "_obj"):
        delattr(node, "_obj")
    with pytest.raises(ImportError):
        _ = node.obj
