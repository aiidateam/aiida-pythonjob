import sys
import types

import cloudpickle
import pytest

from aiida_pythonjob.data.pickled_data import PickledData


def test_roundtrip_basic_types_and_filename():
    data = {"a": 1, "b": [2, 3], "c": ("x", "y")}
    node = PickledData(value=data)

    assert node.get_value() == data
    assert node.value == data

    assert node.get_serialized_value() == cloudpickle.dumps(data)

    assert PickledData.FILENAME == "value.pkl"


def test_roundtrip_function_and_str_representation():
    def multiply(n: int) -> int:
        return n * 2

    node = PickledData(value=multiply)
    fn = node.value
    assert isinstance(fn, types.FunctionType)
    assert fn(3) == 6

    s = str(node)
    assert "PickledData" in s
    assert "multiply" in s  # function name appears


def test_metadata_is_recorded_and_sensible():
    node = PickledData(value=[1, 2, 3])

    attrs = node.base.attributes
    assert attrs.get("python_version") == f"{sys.version_info.major}.{sys.version_info.minor}"
    assert attrs.get("serializer_module") == cloudpickle.__name__
    assert attrs.get("serializer_version") == cloudpickle.__version__
    assert isinstance(attrs.get("pickle_protocol"), int)


def test_setting_value_updates_repository_and_value_property():
    node = PickledData(value="first")
    assert node.value == "first"

    node.value = "second"
    assert node.value == "second"
    assert node.get_serialized_value() == cloudpickle.dumps("second")


def test_error_handling_unpicklingerror(monkeypatch):
    import pickle

    node = PickledData(value="ok")

    def _boom(_bytes):
        raise pickle.UnpicklingError("bad pickle")

    monkeypatch.setattr(cloudpickle, "loads", _boom)

    with pytest.raises(ImportError) as exc:
        _ = node.get_value()
    msg = str(exc.value)
    assert "incompatible pickle protocol" in msg.lower() or "failed to load the pickled value" in msg.lower()


def test_error_handling_valueerror(monkeypatch):
    node = PickledData(value="ok")

    def _boom(_bytes):
        raise ValueError("corrupt stream")

    monkeypatch.setattr(cloudpickle, "loads", _boom)

    with pytest.raises(ImportError) as exc:
        _ = node.get_value()
    assert "failed to load the pickled value" in str(exc.value).lower()


def test_error_handling_missing_module(monkeypatch):
    node = PickledData(value="ok")

    def _boom(_bytes):
        raise ModuleNotFoundError("missing something")

    monkeypatch.setattr(cloudpickle, "loads", _boom)

    with pytest.raises(ImportError) as exc:
        _ = node.get_value()
    msg = str(exc.value).lower()
    assert "missing module" in msg or "failed to load the pickled value" in msg
