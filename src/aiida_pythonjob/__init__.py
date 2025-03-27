"""AiiDA plugin that run Python function on remote computers."""

__version__ = "0.2.4"

from .calculations import PythonJob
from .decorator import pyfunction
from .launch import prepare_pythonjob_inputs
from .parsers import PythonJobParser

__all__ = (
    "PythonJob",
    "pyfunction",
    "PickledData",
    "prepare_pythonjob_inputs",
    "PythonJobParser",
)
