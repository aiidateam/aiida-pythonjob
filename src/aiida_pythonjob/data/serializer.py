from __future__ import annotations

import sys
from importlib.metadata import entry_points
from typing import Any

from aiida import common, orm

from aiida_pythonjob.config import load_config
from aiida_pythonjob.utils import import_from_path

from .deserializer import all_deserializers
from .pickled_data import PickledData


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
    # import time

    # ts = time.time()
    all_serializers = {}
    configs = load_config()
    custom_serializers = configs.get("serializers", {})
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
    data: Any, serializers: dict | None = None, deserializers: dict | None = None, check_value=True
) -> orm.Node:
    """Serialize the data to an AiiDA data node."""
    updated_deserializers = all_deserializers.copy()
    if deserializers is not None:
        updated_deserializers.update(deserializers)

    updated_serializers = all_serializers.copy()
    if serializers is not None:
        updated_serializers.update(serializers)

    if isinstance(data, orm.Data):
        if check_value and not hasattr(data, "value"):
            data_type = type(data)
            ep_key = f"{data_type.__module__}.{data_type.__name__}"
            if ep_key not in updated_deserializers:
                raise ValueError(f"AiiDA data: {ep_key}, does not have a value attribute or deserializer.")
        return data
    elif isinstance(data, common.extendeddicts.AttributeDict):
        # if the data is an AttributeDict, use it directly
        return data
    # if is string with syntax {{}}, this is a port will read data from ctx
    elif isinstance(data, str) and data.startswith("{{") and data.endswith("}}"):
        return data
    # if data is a class instance, get its __module__ and class name as a string
    # for example, an Atoms will have ase.atoms.Atoms
    else:
        data = clean_dict_key(data)
        # try to get the serializer from the entry points
        data_type = type(data)
        ep_key = f"{data_type.__module__}.{data_type.__name__}"
        # search for the key in the entry points
        if ep_key in updated_serializers:
            try:
                serializer = import_from_path(updated_serializers[ep_key])
                new_node = serializer(data)
            except Exception as e:
                raise ValueError(f"Error in serializing {ep_key}: {e}")
            finally:
                # try to save the node to da
                try:
                    new_node.store()
                    return new_node
                except Exception:
                    raise ValueError(f"Error in storing data {ep_key}")
        else:
            # try to serialize the data as a PickledData
            try:
                new_node = PickledData(data)
                new_node.store()
                return new_node
            except Exception as e:
                raise ValueError(f"Error in serializing {ep_key}: {e}")
