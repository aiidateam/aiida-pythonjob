from __future__ import annotations

import sys
import traceback
from importlib.metadata import entry_points
from typing import Any

from aiida import common, orm

from aiida_pythonjob.data.jsonable_data import JsonableData

from .deserializer import all_deserializers
from .utils import import_from_path


def atoms_to_structure_data(structure):
    return orm.StructureData(ase=structure)


def get_serializers_from_entry_points() -> dict:
    # Retrieve the entry points for 'aiida.data' and store them in a dictionary
    eps = entry_points()
    if sys.version_info >= (3, 10):
        group = eps.select(group="aiida.data")
    else:
        group = eps.get("aiida.data", [])
    eps = {}
    for ep in group:
        # split the entry point name by first ".", and check the last part
        key = ep.name.split(".", 1)[-1]
        # skip key without "." because it is not a module name for a data type
        if "." not in key:
            continue
        eps.setdefault(key, [])
        # get the path of the entry point value and replace ":" with "."
        eps[key].append(ep.value.replace(":", "."))
    return eps


def get_serializers() -> dict:
    """Retrieve the serializer from the entry points."""
    from aiida_pythonjob.config import config
    # import time

    # ts = time.time()
    all_serializers = {}
    custom_serializers = config.get("serializers", {})
    eps = get_serializers_from_entry_points()
    # check if there are duplicates
    for key, value in eps.items():
        if len(value) > 1:
            if key not in custom_serializers:
                msg = f"Duplicate entry points for {key}: {value}. You can specify the one to use in the configuration file."  # noqa
                raise ValueError(msg)
        all_serializers[key] = value[0]
    all_serializers.update(custom_serializers)
    # print("Time to get serializer", time.time() - ts)
    return all_serializers


all_serializers = get_serializers()


def serialize_to_aiida_nodes(inputs: dict, serializers: dict | None = None, deserializers: dict | None = None) -> dict:
    """Serialize the inputs to a dictionary of AiiDA data nodes.

    Args:
        inputs (dict): The inputs to be serialized.

    Returns:
        dict: The serialized inputs.
    """
    new_inputs = {}
    # save all kwargs to inputs port
    for key, data in inputs.items():
        new_inputs[key] = general_serializer(data, serializers=serializers, deserializers=deserializers)
    return new_inputs


def clean_dict_key(data):
    """Replace "." with "__dot__" in the keys of a dictionary."""
    if isinstance(data, dict):
        return {k.replace(".", "__dot__"): clean_dict_key(v) for k, v in data.items()}
    return data


def general_serializer(
    data: Any,
    serializers: dict | None = None,
    deserializers: dict | None = None,
    check_value: bool = True,
    store: bool = True,
) -> orm.Node:
    """
    Attempt to serialize the data to an AiiDA data node based on the preference from `config`:
      1) AiiDA data only, 2) JSON-serializable, 3) fallback to PickledData (if allowed).
    """
    from aiida_pythonjob.config import config

    # Merge user-provided config with defaults
    allow_json = config.get("allow_json", True)
    allow_pickle = config.get("allow_pickle", False)

    updated_deserializers = all_deserializers.copy()
    if deserializers is not None:
        updated_deserializers.update(deserializers)

    updated_serializers = all_serializers.copy()
    if serializers is not None:
        updated_serializers.update(serializers)

    # 1) If it is already an AiiDA node, just return it
    if isinstance(data, orm.Data):
        if check_value and not hasattr(data, "value"):
            data_type = type(data)
            ep_key = f"{data_type.__module__}.{data_type.__name__}"
            if ep_key not in updated_deserializers:
                raise ValueError(f"AiiDA data: {ep_key}, does not have a `value` attribute or deserializer.")
        return data
    elif isinstance(data, common.extendeddicts.AttributeDict):
        # if the data is an AttributeDict, use it directly
        return data

    # 3) check entry point
    data_type = type(data)
    ep_key = f"{data_type.__module__}.{data_type.__name__}"
    if ep_key in updated_serializers:
        try:
            serializer = import_from_path(updated_serializers[ep_key])
            new_node = serializer(data)
            if store:
                new_node.store()
            return new_node
        except Exception:
            error_traceback = traceback.format_exc()
            raise ValueError(f"Error in serializing {ep_key}: {error_traceback}")

    #    check if we can JSON-serialize the data
    if allow_json:
        try:
            node = JsonableData(data)
            if store:
                node.store()
            return node
        except (TypeError, ValueError):
            print(f"Error in JSON-serializing {type(data).__name__}")

    # fallback to pickling
    if allow_pickle:
        from .pickled_data import PickledData

        try:
            new_node = PickledData(data)
            if store:
                new_node.store()
            return new_node
        except Exception as e:
            raise ValueError(f"Error in pickling {type(data).__name__}: {e}")

    raise ValueError(
        f"Cannot serialize type={type(data).__name__}. No suitable method found "
        f"(json_allowed={allow_json}, pickle_allowed={allow_pickle})."
    )
