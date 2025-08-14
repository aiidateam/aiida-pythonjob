from __future__ import annotations

from typing import Any, Dict, Optional


def _is_namespace_identifier(identifier: str, *, ns_id: str = "node_graph.namespace") -> bool:
    return identifier == ns_id or identifier.endswith(".namespace")


def _socket_meta_required(sock: Dict[str, Any]) -> Optional[bool]:
    meta = sock.get("metadata", {}) or {}
    if "required" in meta and meta["required"] is not None:
        return bool(meta["required"])
    if sock.get("property", {}).get("default", None) is not None:
        return False
    return None


def _socket_meta_help(sock: Dict[str, Any]) -> Optional[str]:
    return (sock.get("metadata") or {}).get("help")


def _socket_to_port_object(sock: Dict[str, Any]) -> Dict[str, Any]:
    ident = sock.get("identifier", "")
    meta = sock.get("metadata", {}) or {}

    if _is_namespace_identifier(ident):
        ports: Dict[str, Dict[str, Any]] = {}
        for name, child in (sock.get("sockets") or {}).items():
            ports[name] = _socket_to_port_object(child)

        obj: Dict[str, Any] = {"identifier": "NAMESPACE", "ports": ports}

        if meta.get("dynamic"):
            obj["dynamic"] = True
            # Prefer full nested item namespace when present
            if isinstance(meta.get("item"), dict):
                obj["item"] = _socket_to_port_object(meta["item"])
            else:
                # fallback: identifier only (treat namespace id as empty namespace)
                item_ident = meta.get("item_identifier", "node_graph.any")
                obj["item"] = (
                    {"identifier": "NAMESPACE", "ports": {}}
                    if _is_namespace_identifier(item_ident)
                    else {"identifier": "ANY"}
                )

        h = _socket_meta_help(sock)
        if h:
            obj["help"] = h
        req = _socket_meta_required(sock)
        if req is not None:
            obj["required"] = req
        return obj

    # Leaf
    obj = {"identifier": "ANY"}
    h = _socket_meta_help(sock)
    if h:
        obj["help"] = h
    req = _socket_meta_required(sock)
    if req is not None:
        obj["required"] = req
    return obj


def _port_obj_to_list_shape(obj: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively convert {'ports': {name: obj}} → {'ports': [{'name': name, ...}] }."""
    if obj.get("identifier") == "NAMESPACE":
        mapped = []
        for name, child in (obj.get("ports") or {}).items():
            mapped.append({"name": name, **_port_obj_to_list_shape(child)})
        out = {"identifier": "NAMESPACE", "ports": mapped}
        if obj.get("dynamic"):
            out["dynamic"] = True
        if "item" in obj:
            out["item"] = _port_obj_to_list_shape(obj["item"])
        # carry optional fields
        if "help" in obj:
            out["help"] = obj["help"]
        if "required" in obj:
            out["required"] = obj["required"]
        return out
    # leaf: nothing to convert
    return {k: v for k, v in obj.items() if k != "ports"}


def inputs_sockets_to_ports(node_inputs: Dict[str, Any], *, as_list: bool = False) -> Dict[str, Any]:
    ports_map: Dict[str, Dict[str, Any]] = {}
    for name, sock in (node_inputs.get("sockets") or {}).items():
        ports_map[name] = _socket_to_port_object(sock)

    schema: Dict[str, Any] = {"name": "inputs", "identifier": "NAMESPACE", "ports": ports_map}
    meta = node_inputs.get("metadata", {}) or {}
    if meta.get("dynamic"):
        schema["dynamic"] = True
        if isinstance(meta.get("item"), dict):
            schema["item"] = _socket_to_port_object(meta["item"])
        else:
            item_ident = meta.get("item_identifier", "node_graph.any")
            schema["item"] = (
                {"identifier": "NAMESPACE", "ports": {}}
                if _is_namespace_identifier(item_ident)
                else {"identifier": "ANY"}
            )

    if as_list:
        schema["ports"] = [{"name": n, **_port_obj_to_list_shape(p)} for n, p in ports_map.items()]
        if "item" in schema:
            schema["item"] = _port_obj_to_list_shape(schema["item"])
    return schema


def outputs_sockets_to_ports(node_outputs: Dict[str, Any], *, as_list: bool = False) -> Dict[str, Any]:
    meta = node_outputs.get("metadata", {}) or {}
    sockets = node_outputs.get("sockets") or {}

    ports_map = {name: _socket_to_port_object(sock) for name, sock in sockets.items()}
    schema: Dict[str, Any] = {"name": "outputs", "identifier": "NAMESPACE", "ports": ports_map}
    if meta.get("dynamic"):
        schema["dynamic"] = True
        if isinstance(meta.get("item"), dict):
            schema["item"] = _socket_to_port_object(meta["item"])
        else:
            item_ident = meta.get("item_identifier", "node_graph.any")
            schema["item"] = (
                {"identifier": "NAMESPACE", "ports": {}}
                if _is_namespace_identifier(item_ident)
                else {"identifier": "ANY"}
            )
    if as_list:
        schema["ports"] = [{"name": n, **_port_obj_to_list_shape(p)} for n, p in ports_map.items()]
        if "item" in schema:
            schema["item"] = _port_obj_to_list_shape(schema["item"])
    return schema
