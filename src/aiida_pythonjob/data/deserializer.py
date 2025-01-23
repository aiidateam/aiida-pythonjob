from typing import Any

from aiida import common, orm

from aiida_pythonjob.config import load_config

builtin_deserializers = {
    "aiida.orm.nodes.data.list.List": "aiida_pythonjob.data.deserializer.list_data_to_list",
    "aiida.orm.nodes.data.dict.Dict": "aiida_pythonjob.data.deserializer.dict_data_to_dict",
}


def list_data_to_list(data):
    return data.get_list()


def dict_data_to_dict(data):
    return data.get_dict()


def structure_data_to_atoms(structure):
    return structure.get_ase()


def get_deserializer() -> dict:
    """Retrieve the serializer from the entry points."""
    configs = load_config()
    custom_deserializers = configs.get("deserializers", {})
    deserializers = builtin_deserializers.copy()
    deserializers.update(custom_deserializers)
    return deserializers


eps_deserializers = get_deserializer()


def deserialize_to_raw_python_data(datas: dict) -> dict:
    """Deserialize the datas to a dictionary of raw Python data.

    Args:
        datas (dict): The datas to be deserialized.

    Returns:
        dict: The deserialized datas.
    """
    new_datas = {}
    # save all kwargs to inputs port
    for key, data in datas.items():
        new_datas[key] = general_deserializer(data)
    return new_datas


def general_deserializer(data: Any) -> orm.Node:
    """Deserialize the AiiDA data node to an raw Python data."""
    import importlib

    if isinstance(data, orm.Data):
        if hasattr(data, "value"):
            return getattr(data, "value")
        data_type = type(data)
        ep_key = f"{data_type.__module__}.{data_type.__name__}"
        if ep_key in eps_deserializers:
            module_name, deserializer_name = eps_deserializers[ep_key].rsplit(".", 1)
            module = importlib.import_module(module_name)
            deserializer = getattr(module, deserializer_name)
            return deserializer(data)
        else:
            raise ValueError(f"AiiDA data: {ep_key}, does not have a value attribute or deserializer.")
    elif isinstance(data, (common.extendeddicts.AttributeDict, dict)):
        # if the data is an AttributeDict, use it directly
        return {k: general_deserializer(v) for k, v in data.items()}
