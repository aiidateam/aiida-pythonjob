"""AiiDA plugin that run Python function on remote computers."""

__version__ = "0.1.4"

from .calculations import PythonJob
from .launch import prepare_pythonjob_inputs
from .parsers import PythonJobParser

__all__ = (
    "PythonJob",
    "PickledData",
    "prepare_pythonjob_inputs",
    "PythonJobParser",
)
