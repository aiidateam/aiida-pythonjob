"""Proess to run a Python function locally"""

from __future__ import annotations

import inspect
import typing as t

from aiida.common.lang import override
from aiida.engine import Process, ProcessSpec
from aiida.engine.processes.exit_code import ExitCode
from aiida.orm import (
    CalcFunctionNode,
    Data,
    Dict,
    List,
    Str,
    to_aiida_type,
)

__all__ = ("PyFunction",)


# The following code is modified from the aiida-core.engine.processes.functions module


# TODO: because aiida-core cmds hardcode `CalcFunctionNode`, I hardcoded `CalcFunctionNode` here.
# In principle, aiida-core should support subclassing
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
        if self._func is None:
            self._func = self.inputs.function_data.pickled_function.value
        return self._func

    @classmethod
    def define(cls, spec: ProcessSpec) -> None:  # type: ignore[override]
        """Define the process specification, including its inputs, outputs and known exit codes."""
        super().define(spec)
        spec.input_namespace("function_data")
        spec.input("function_data.name", valid_type=Str, serializer=to_aiida_type)
        spec.input("function_data.source_code", valid_type=Str, serializer=to_aiida_type, required=False)
        spec.input("function_data.outputs", valid_type=List, serializer=to_aiida_type, required=False)
        spec.input("function_data.pickled_function", valid_type=Data, required=False)
        spec.input("function_data.mode", valid_type=Str, serializer=to_aiida_type, required=False)
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

    def get_function_name(self) -> str:
        """Return the name of the function to run."""
        if "name" in self.inputs.function_data:
            name = self.inputs.function_data.name.value
        else:
            try:
                name = self.inputs.function_data.pickled_function.value.__name__
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

        from aiida_pythonjob.data.deserializer import deserialize_to_raw_python_data

        # The following conditional is required for the caching to properly work. Even if the source node has a process
        # state of `Finished` the cached process will still enter the running state. The process state will have then
        # been overridden by the engine to `Running` so we cannot check that, but if the `exit_status` is anything other
        # than `None`, it should mean this node was taken from the cache, so the process should not be rerun.
        if self.node.exit_status is not None:
            return ExitCode(self.node.exit_status, self.node.exit_message)

        # Now the original functions arguments need to be reconstructed from the inputs to the process, as they were
        # passed to the original function call. To do so, all positional parameters are popped from the inputs
        # dictionary and added to the positional arguments list.
        args = []
        kwargs: dict[str, Data] = {}
        inputs = dict(self.inputs.function_inputs or {})
        inputs.pop("serializers", None)
        inputs.pop("deserializers", None)
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

        for name, parameter in inspect.signature(self.func).parameters.items():
            if parameter.kind in [parameter.POSITIONAL_ONLY, parameter.POSITIONAL_OR_KEYWORD]:
                if name in inputs:
                    args.append(inputs.pop(name))
            elif parameter.kind is parameter.VAR_POSITIONAL:
                for key in [key for key in inputs.keys() if key.startswith(f"{name}_")]:
                    args.append(inputs.pop(key))

        # Any inputs that correspond to metadata ports were not part of the original function signature but were added
        # by the process function decorator, so these have to be removed.
        for key in [key for key, port in self.spec().inputs.items() if port.is_metadata]:
            inputs.pop(key, None)

        # The remaining inputs have to be keyword arguments.
        kwargs.update(**inputs)
        raw_args = [deserialize_to_raw_python_data(x, deserializers=self.deserializers) for x in args]
        raw_kwargs = {k: deserialize_to_raw_python_data(v, deserializers=self.deserializers) for k, v in kwargs.items()}

        results = self.func(*raw_args, **raw_kwargs)

        # Read function_outputs specification
        if "outputs" in self.inputs.function_data:
            function_outputs = self.node.inputs.function_data.outputs.get_list()
        else:
            function_outputs = [{"name": "result"}]
        self.output_list = function_outputs

        # If nested outputs like "add_multiply.add", keep only top-level
        top_level_output_list = [output for output in self.output_list if "." not in output["name"]]
        if isinstance(results, tuple):
            if len(top_level_output_list) != len(results):
                return self.exit_codes.ERROR_RESULT_OUTPUT_MISMATCH
            for i in range(len(top_level_output_list)):
                top_level_output_list[i]["value"] = self.serialize_output(results[i], top_level_output_list[i])
        elif isinstance(results, dict):
            # pop the exit code if it exists inside the dictionary
            exit_code = results.pop("exit_code", None)
            if exit_code:
                # If there's an exit_code, handle it (dict or int)
                if isinstance(exit_code, dict):
                    exit_code = ExitCode(exit_code["status"], exit_code["message"])
                elif isinstance(exit_code, int):
                    exit_code = ExitCode(exit_code)
                if exit_code.status != 0:
                    return exit_code
            if len(top_level_output_list) == 1:
                # If output name in results, use it
                if top_level_output_list[0]["name"] in results:
                    top_level_output_list[0]["value"] = self.serialize_output(
                        results.pop(top_level_output_list[0]["name"]),
                        top_level_output_list[0],
                    )
                    # If there are any extra keys in `results`, log a warning
                    if len(results) > 0:
                        self.logger.warning(
                            f"Found extra results that are not included in the output: {results.keys()}"
                        )
                else:
                    # Otherwise assume the entire dict is the single output
                    top_level_output_list[0]["value"] = self.serialize_output(results, top_level_output_list[0])
            elif len(top_level_output_list) > 1:
                # Match each top-level output by name
                for output in top_level_output_list:
                    if output["name"] not in results:
                        if output.get("required", True):
                            return self.exit_codes.ERROR_MISSING_OUTPUT
                    else:
                        output["value"] = self.serialize_output(results.pop(output["name"]), output)
                # Any remaining results are unaccounted for -> log a warning
                if len(results) > 0:
                    self.logger.warning(f"Found extra results that are not included in the output: {results.keys()}")
        elif len(top_level_output_list) == 1:
            # Single top-level output, single result
            top_level_output_list[0]["value"] = self.serialize_output(results, top_level_output_list[0])
        else:
            return self.exit_codes.ERROR_RESULT_OUTPUT_MISMATCH
        # Store the outputs
        for output in top_level_output_list:
            self.out(output["name"], output["value"])

        return ExitCode()

    def find_output(self, name):
        """Find the output spec with the given name."""
        for output in self.output_list:
            if output["name"] == name:
                return output
        return None

    def serialize_output(self, result, output):
        """Serialize outputs."""
        from aiida_pythonjob.data.serializer import general_serializer

        name = output["name"]
        if output.get("identifier", "Any").upper() == "NAMESPACE":
            if isinstance(result, dict):
                serialized_result = {}
                for key, value in result.items():
                    full_name = f"{name}.{key}"
                    full_name_output = self.find_output(full_name)
                    if full_name_output and full_name_output.get("identifier", "Any").upper() == "NAMESPACE":
                        serialized_result[key] = self.serialize_output(value, full_name_output)
                    else:
                        serialized_result[key] = general_serializer(value, serializers=self.serializers, store=False)
                return serialized_result
            else:
                self.logger.error(f"Expected a dict for namespace '{name}', got {type(result)}.")
                return self.exit_codes.ERROR_INVALID_OUTPUT
        else:
            return general_serializer(result, serializers=self.serializers, store=False)
