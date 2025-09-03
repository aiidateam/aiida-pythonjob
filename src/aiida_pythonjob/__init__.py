"""AiiDA plugin that run Python function on remote computers."""

__version__ = "0.3.2"

from node_graph import socket_spec as spec

from .calculations import PyFunction, PythonJob
from .decorator import pyfunction
from .launch import prepare_pythonjob_inputs
from .parsers import PythonJobParser

__all__ = (
    "PythonJob",
    "PyFunction",
    "pyfunction",
    "PickledData",
    "prepare_pythonjob_inputs",
    "PythonJobParser",
    "spec",
)
