"""Proess to run a Python function locally"""

from __future__ import annotations

import traceback
import typing as t

import cloudpickle
import plumpy
from aiida.common.lang import override
from aiida.engine import Process, ProcessSpec
from aiida.engine.processes.exit_code import ExitCode
from aiida.orm import CalcFunctionNode
from node_graph.socket_spec import SocketSpec

from aiida_pythonjob.calculations.common import (
    ATTR_DESERIALIZERS,
    ATTR_OUTPUTS_SPEC,
    ATTR_SERIALIZERS,
    FunctionProcessMixin,
    add_common_function_io,
)
from aiida_pythonjob.data.deserializer import deserialize_to_raw_python_data
from aiida_pythonjob.parsers.utils import parse_outputs

__all__ = ("PyFunction",)


class PyFunction(FunctionProcessMixin, Process):
    """Run a Python function in-process, using :class:`SocketSpec` for I/O."""

    _node_class = CalcFunctionNode
    label_template = "{name}"
    default_name = "anonymous_function"

    def __init__(self, *args, **kwargs) -> None:
        if kwargs.get("enable_persistence", False):
            raise RuntimeError("Cannot persist a function process")
        super().__init__(enable_persistence=False, *args, **kwargs)  # type: ignore[misc]
        self._func = None

    @override
    def load_instance_state(
        self, saved_state: t.MutableMapping[str, t.Any], load_context: plumpy.persistence.LoadSaveContext
    ) -> None:
        """Load the instance state (restore pickled function)."""
        super().load_instance_state(saved_state, load_context)
        self._func = cloudpickle.loads(self.inputs.function_data.pickled_function)

    @property
    def func(self) -> t.Callable[..., t.Any]:
        if self._func is None:
            self._func = cloudpickle.loads(self.inputs.function_data.pickled_function)
        return self._func

    def _extract_declared_name(self) -> str | None:
        name = super()._extract_declared_name()
        if name:
            return name
        try:
            return self.func.__name__
        except Exception:
            return None

    @classmethod
    def define(cls, spec: ProcessSpec) -> None:  # type: ignore[override]
        """Define inputs/outputs and exit codes."""
        super().define(spec)
        add_common_function_io(spec)
        spec.inputs.dynamic = True
        spec.outputs.dynamic = True
        spec.exit_code(
            323,
            "ERROR_FUNCTION_EXECUTION_FAILED",
            invalidates_cache=True,
            message="Function execution failed.\n{exception}\n{traceback}",
        )

    @override
    def _setup_db_record(self) -> None:
        super()._setup_db_record()
        self.node.store_source_info(self.func)

    def execute(self) -> dict[str, t.Any] | None:
        """Mirror calcfunction behavior: unwrap single-output dicts to a bare value."""
        result = super().execute()
        if result and len(result) == 1 and self.SINGLE_OUTPUT_LINKNAME in result:
            return result[self.SINGLE_OUTPUT_LINKNAME]
        return result

    @override
    def run(self) -> ExitCode | None:
        # Respect caching semantics (from aiida-core calcfunction implementation)
        if self.node.exit_status is not None:
            return ExitCode(self.node.exit_status, self.node.exit_message)

        # Deserialize inputs
        try:
            inputs = dict(self.inputs.function_inputs or {})
            deserializers = self.node.base.attributes.get(ATTR_DESERIALIZERS, {})
            inputs = deserialize_to_raw_python_data(inputs, deserializers=deserializers)
        except Exception as exception:
            return self.exit_codes.ERROR_DESERIALIZE_INPUTS_FAILED.format(
                exception=str(exception), traceback=traceback.format_exc()
            )

        # Execute function
        try:
            results = self.func(**inputs)
        except Exception as exception:
            return self.exit_codes.ERROR_FUNCTION_EXECUTION_FAILED.format(
                exception=str(exception), traceback=traceback.format_exc()
            )

        # Parse & attach outputs
        outputs_spec = SocketSpec.from_dict(self.node.base.attributes.get(ATTR_OUTPUTS_SPEC) or {})
        serializers = self.node.base.attributes.get(ATTR_SERIALIZERS, {})
        outputs, exit_code = parse_outputs(
            results,
            output_spec=outputs_spec,
            exit_codes=self.exit_codes,
            logger=self.logger,
            serializers=serializers,
        )
        if exit_code:
            return exit_code

        for name, value in (outputs or {}).items():
            self.out(name, value)
        return ExitCode()
