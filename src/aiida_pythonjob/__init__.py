"""AiiDA plugin that run Python function on remote computers."""

__version__ = "0.3.5"

from node_graph import socket_spec as spec

from .calculations import PyFunction, PythonJob
from .decorator import pyfunction
from .launch import prepare_pyfunction_inputs, prepare_pythonjob_inputs
from .parsers import PythonJobParser

__all__ = (
    "PickledData",
    "PyFunction",
    "PythonJob",
    "PythonJobParser",
    "prepare_pyfunction_inputs",
    "prepare_pythonjob_inputs",
    "pyfunction",
    "spec",
)
