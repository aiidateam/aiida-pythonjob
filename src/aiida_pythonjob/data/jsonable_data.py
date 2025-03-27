import importlib
import json
import typing

import numpy as np
from aiida import orm

__all__ = ("JsonableData",)


class JsonableData(orm.Data):
    """
    A more flexible version of `JsonableData` that:
      1) Detects multiple possible "to-dict" methods (as_dict, to_dict, todict, etc.)
      2) Converts non-JSON-serializable items (e.g. NumPy arrays) into JSON-friendly formats
      3) Rebuilds the object via one of various possible "from-dict" methods.
      4) Safely handles `None`.
    """

    _DICT_METHODS = ["as_dict", "to_dict", "todict", "asdict"]

    _FROM_DICT_METHODS = ["from_dict", "fromdict"]

    def __init__(self, obj: typing.Any, *args, **kwargs):
        """
        Construct the node for the to-be-wrapped object.
        If `obj` is None, we simply store a flag indicating that.
        """
        super().__init__(*args, **kwargs)

        if obj is None:
            self.base.attributes.set("is_none", True)
            self._obj = None
            return

        dictionary = self._extract_dict(obj)

        dictionary.setdefault("@class", obj.__class__.__name__)
        dictionary.setdefault("@module", obj.__class__.__module__)

        dictionary = self._make_jsonable(dictionary)

        # Because some float constants (NaN, inf) can cause trouble with PostgreSQL JSON,
        # we do a round-trip through Python's JSON module with a custom parse_constant
        try:
            serialized = json.loads(json.dumps(dictionary), parse_constant=lambda x: x)
        except TypeError as exc:
            raise TypeError(f"Object `{obj}` cannot be fully JSON-serialized.") from exc

        self.base.attributes.set_many(serialized)
        self._obj = obj

    def _extract_dict(self, obj: typing.Any) -> dict:
        """
        Attempt to call one of the recognized "to-dict" style methods on `obj` in sequence.
        """
        for method_name in self._DICT_METHODS:
            method = getattr(obj, method_name, None)
            if callable(method):
                return method()

        raise TypeError(
            f"Object `{obj}` does not have any of the following dictionary-conversion " f"methods: {self._DICT_METHODS}"
        )

    def _make_jsonable(self, data: typing.Any) -> typing.Any:
        """
        Recursively walk `data`. Convert anything that is not JSON-serializable
        into JSON-friendly structures (e.g. convert NumPy arrays to lists).
        """
        if isinstance(data, dict):
            return {k: self._make_jsonable(v) for k, v in data.items()}

        elif isinstance(data, list):
            return [self._make_jsonable(v) for v in data]

        elif isinstance(data, tuple):
            return tuple(self._make_jsonable(v) for v in data)

        elif isinstance(data, np.ndarray):
            return data.tolist()

        elif isinstance(data, (np.generic,)):  # e.g. np.int64, np.float32
            return data.item()

        return data

    @classmethod
    def _deserialize_float_constants(cls, data: typing.Any) -> typing.Any:
        """
        Handle Infinity, -Infinity, NaN from round-tripped JSON.
        """
        if isinstance(data, dict):
            return {k: cls._deserialize_float_constants(v) for k, v in data.items()}
        if isinstance(data, list):
            return [cls._deserialize_float_constants(v) for v in data]
        if data == "Infinity":
            return float("inf")
        if data == "-Infinity":
            return -float("inf")
        if data == "NaN":
            return float("nan")
        return data

    def _get_object(self) -> typing.Any:
        """
        Return the wrapped Python object. If not cached in `_obj`,
        we reconstruct it by calling one of the recognized 'from-dict' methods.
        """
        if self.base.attributes.get("is_none", False):
            self._obj = None
            return None

        if hasattr(self, "_obj"):
            return self._obj

        attributes = self.base.attributes.all
        class_name = attributes.pop("@class")
        module_name = attributes.pop("@module")

        try:
            module = importlib.import_module(module_name)
        except ImportError as exc:
            raise ImportError(f"The object's module `{module_name}` cannot be imported.") from exc

        try:
            cls_ = getattr(module, class_name)
        except AttributeError as exc:
            raise ImportError(f"The module `{module_name}` does not contain the class `{class_name}`.") from exc

        attributes["@class"] = class_name
        attributes["@module"] = module_name

        deserialized = self._deserialize_float_constants(attributes)
        self._obj = self._rebuild_object(cls_, deserialized)
        return self._obj

    def _rebuild_object(self, cls_: typing.Any, attributes: dict) -> typing.Any:
        """
        Attempt to reconstruct an object of type `cls_` from `attributes`.
        """
        for method_name in self._FROM_DICT_METHODS:
            fromdict_method = getattr(cls_, method_name, None)
            if callable(fromdict_method):
                return fromdict_method(attributes)

        try:
            return cls_(**attributes)
        except TypeError:
            raise TypeError(
                f"Cannot rebuild object of type `{cls_}`: no suitable from-dict method "
                f"({self._FROM_DICT_METHODS}) nor constructor that accepts these attributes."
            )

    @property
    def obj(self) -> typing.Any:
        """Return the wrapped Python object (or None if `is_none` is set)."""
        return self._get_object()

    @property
    def value(self) -> typing.Any:
        """Alias for `.obj`."""
        return self._get_object()
