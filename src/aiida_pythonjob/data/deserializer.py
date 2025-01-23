from __future__ import annotations

from typing import Any

from aiida import common, orm

from aiida_pythonjob.config import load_config

builtin_deserializers = {
    "aiida.orm.nodes.data.list.List": "aiida_pythonjob.data.deserializer.list_data_to_list",
    "aiida.orm.nodes.data.dict.Dict": "aiida_pythonjob.data.deserializer.dict_data_to_dict",
    "aiida.orm.nodes.data.structure.StructureData": "aiida_pythonjob.data.deserializer.structure_data_to_atoms",
}


def generate_aiida_node_deserializer(data: orm.Node) -> dict:
    if isinstance(data, orm.Data):
        return data.backend_entity.attributes
    elif isinstance(data, (common.extendeddicts.AttributeDict, dict)):
        # if the data is an AttributeDict, use it directly
        return {k: generate_aiida_node_deserializer(v) for k, v in data.items()}


def list_data_to_list(data):
    return data.get_list()


def dict_data_to_dict(data):
    return data.get_dict()


def structure_data_to_atoms(structure):
    return structure.get_ase()


def structure_data_to_pymatgen(structure):
    return structure.get_pymatgen()


def get_deserializer() -> dict:
    """Retrieve the serializer from the entry points."""
    configs = load_config()
    custom_deserializers = configs.get("deserializers", {})
    deserializers = builtin_deserializers.copy()
    deserializers.update(custom_deserializers)
    return deserializers


eps_deserializers = get_deserializer()


def deserialize_to_raw_python_data(data: orm.Node, deserializers: dict | None = None) -> Any:
    """Deserialize the AiiDA data node to an raw Python data."""
    import importlib

    all_deserializers = eps_deserializers.copy()

    if deserializers is not None:
        all_deserializers.update(deserializers)

    if isinstance(data, orm.Data):
        if hasattr(data, "value"):
            return getattr(data, "value")
        data_type = type(data)
        ep_key = f"{data_type.__module__}.{data_type.__name__}"
        if ep_key in all_deserializers:
            module_name, deserializer_name = all_deserializers[ep_key].rsplit(".", 1)
            module = importlib.import_module(module_name)
            deserializer = getattr(module, deserializer_name)
            return deserializer(data)
        else:
            raise ValueError(f"AiiDA data: {ep_key}, does not have a value attribute or deserializer.")
    elif isinstance(data, (common.extendeddicts.AttributeDict, dict)):
        # if the data is an AttributeDict, use it directly
        return {k: deserialize_to_raw_python_data(v, deserializers=deserializers) for k, v in data.items()}
