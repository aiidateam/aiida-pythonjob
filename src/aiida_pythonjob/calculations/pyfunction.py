"""Proess to run a Python function locally"""

from __future__ import annotations

import traceback
import typing as t

from aiida.common.lang import override
from aiida.engine import Process, ProcessSpec
from aiida.engine.processes.exit_code import ExitCode
from aiida.orm import (
    CalcFunctionNode,
    Data,
    Dict,
    Str,
    to_aiida_type,
)

__all__ = ("PyFunction",)


class PyFunction(Process):
    """"""

    _node_class = CalcFunctionNode

    def __init__(self, *args, **kwargs) -> None:
        if kwargs.get("enable_persistence", False):
            raise RuntimeError("Cannot persist a function process")
        super().__init__(enable_persistence=False, *args, **kwargs)  # type: ignore[misc]
        self._func = None

    @property
    def func(self) -> t.Callable[..., t.Any]:
        import cloudpickle

        if self._func is None:
            self._func = cloudpickle.loads(self.inputs.function_data.pickled_function)
        return self._func

    @classmethod
    def define(cls, spec: ProcessSpec) -> None:  # type: ignore[override]
        """Define the process specification, including its inputs, outputs and known exit codes."""
        super().define(spec)
        spec.input_namespace("function_data", dynamic=True, required=True)
        spec.input("function_data.output_ports", valid_type=Dict, serializer=to_aiida_type, required=False)
        spec.input("function_data.input_ports", valid_type=Dict, serializer=to_aiida_type, required=False)
        spec.input("process_label", valid_type=Str, serializer=to_aiida_type, required=False)
        spec.input_namespace("function_inputs", valid_type=Data, required=False)
        spec.input(
            "deserializers",
            valid_type=Dict,
            default=None,
            required=False,
            serializer=to_aiida_type,
            help="The deserializers to convert the input AiiDA data nodes to raw Python data.",
        )
        spec.input(
            "serializers",
            valid_type=Dict,
            default=None,
            required=False,
            serializer=to_aiida_type,
            help="The serializers to convert the raw Python data to AiiDA data nodes.",
        )
        spec.inputs.dynamic = True
        spec.outputs.dynamic = True
        spec.exit_code(
            320,
            "ERROR_INVALID_OUTPUT",
            invalidates_cache=True,
            message="The output file contains invalid output.",
        )
        spec.exit_code(
            321,
            "ERROR_RESULT_OUTPUT_MISMATCH",
            invalidates_cache=True,
            message="The number of results does not match the number of outputs.",
        )
        spec.exit_code(
            323,
            "ERROR_DESERIALIZE_INPUTS_FAILED",
            invalidates_cache=True,
            message="Failed to unpickle inputs.\n{exception}\n{traceback}",
        )
        spec.exit_code(
            325,
            "ERROR_FUNCTION_EXECUTION_FAILED",
            invalidates_cache=True,
            message="Function execution failed.\n{exception}\n{traceback}",
        )

    def get_function_name(self) -> str:
        """Return the name of the function to run."""
        if "name" in self.inputs.function_data:
            name = self.inputs.function_data.name
        else:
            try:
                name = self.func.__name__
            except AttributeError:
                # If a user doesn't specify name, fallback to something generic
                name = "anonymous_function"
        return name

    def _build_process_label(self) -> str:
        """Use the function name or an explicit label as the process label."""
        if "process_label" in self.inputs:
            return self.inputs.process_label.value
        else:
            name = self.get_function_name()
            return f"{name}"

    @override
    def _setup_db_record(self) -> None:
        """Set up the database record for the process."""
        super()._setup_db_record()
        self.node.store_source_info(self.func)

    def execute(self) -> dict[str, t.Any] | None:
        """Execute the process."""
        result = super().execute()

        # FunctionProcesses can return a single value as output, and not a dictionary, so we should also return that
        if result and len(result) == 1 and self.SINGLE_OUTPUT_LINKNAME in result:
            return result[self.SINGLE_OUTPUT_LINKNAME]

        return result

    @override
    def run(self) -> ExitCode | None:
        """Run the process."""
        from aiida_pythonjob.utils import deserialize_ports

        # The following conditional is required for the caching to properly work.
        # From the the calcfunction implementation in aiida-core
        if self.node.exit_status is not None:
            return ExitCode(self.node.exit_status, self.node.exit_message)

        # load custom serializers
        if "serializers" in self.node.inputs and self.node.inputs.serializers:
            serializers = self.node.inputs.serializers.get_dict()
            # replace "__dot__" with "." in the keys
            self.serializers = {k.replace("__dot__", "."): v for k, v in serializers.items()}
        else:
            self.serializers = None
        if "deserializers" in self.node.inputs and self.node.inputs.deserializers:
            deserializers = self.node.inputs.deserializers.get_dict()
            # replace "__dot__" with "." in the keys
            self.deserializers = {k.replace("__dot__", "."): v for k, v in deserializers.items()}
        else:
            self.deserializers = None

        inputs = dict(self.inputs.function_inputs or {})
        try:
            inputs = deserialize_ports(
                serialized_data=inputs,
                port_schema=self.node.inputs.function_data.input_ports.get_dict(),
                deserializers=self.deserializers,
            )
        except Exception as exception:
            exception_message = str(exception)
            traceback_str = traceback.format_exc()
            return self.exit_codes.ERROR_DESERIALIZE_INPUTS_FAILED.format(
                exception=exception_message, traceback=traceback_str
            )

        try:
            results = self.func(**inputs)
        except Exception as exception:
            exception_message = str(exception)
            traceback_str = traceback.format_exc()
            return self.exit_codes.ERROR_FUNCTION_EXECUTION_FAILED.format(
                exception=exception_message, traceback=traceback_str
            )
        return self.parse(results)

    def parse(self, results):
        """Parse the results of the function."""
        from aiida_pythonjob.utils import parse_outputs

        self.output_ports = self.node.inputs.function_data.output_ports.get_dict()
        exit_code = parse_outputs(
            results,
            output_ports=self.output_ports,
            exit_codes=self.exit_codes,
            logger=self.logger,
            serializers=self.serializers,
        )
        if exit_code:
            return exit_code
        # Store the outputs
        for output in self.output_ports["ports"]:
            self.out(output["name"], output["value"])

        return ExitCode()
