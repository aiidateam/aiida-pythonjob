"""Parser for an `PythonJob` job."""

import json

from aiida.engine import ExitCode
from aiida.parsers.parser import Parser

# Map error_type from script.py to exit code label
ERROR_TYPE_TO_EXIT_CODE = {
    "IMPORT_CLOUDPICKLE_FAILED": "ERROR_IMPORT_CLOUDPICKLE_FAILED",
    "UNPICKLE_INPUTS_FAILED": "ERROR_UNPICKLE_INPUTS_FAILED",
    "UNPICKLE_FUNCTION_FAILED": "ERROR_UNPICKLE_FUNCTION_FAILED",
    "FUNCTION_EXECUTION_FAILED": "ERROR_FUNCTION_EXECUTION_FAILED",
    "PICKLE_RESULTS_FAILED": "ERROR_PICKLE_RESULTS_FAILED",
}


class PythonJobParser(Parser):
    """Parser for an `PythonJob` job."""

    def parse(self, **kwargs):
        import pickle

        # Read function_outputs specification
        if "outputs" in self.node.inputs.function_data:
            function_outputs = self.node.inputs.function_data.outputs.get_list()
        else:
            function_outputs = [{"name": "result"}]
        self.output_list = function_outputs

        # load custom serializers
        if "serializers" in self.node.inputs and self.node.inputs.serializers:
            serializers = self.node.inputs.serializers.get_dict()
            # replace "__dot__" with "." in the keys
            self.serializers = {k.replace("__dot__", "."): v for k, v in serializers.items()}
        else:
            self.serializers = None

        # If nested outputs like "add_multiply.add", keep only top-level
        top_level_output_list = [output for output in self.output_list if "." not in output["name"]]

        # 1) Read _error.json
        error_data = {}
        try:
            with self.retrieved.base.repository.open("_error.json", "r") as ef:
                error_data = json.load(ef)
        except OSError:
            # No _error.json file found
            pass
        except json.JSONDecodeError as exc:
            self.logger.error(f"Error reading _error.json: {exc}")
            return self.exit_codes.ERROR_INVALID_OUTPUT  # or a different exit code

        # If error_data is non-empty, we have an error from the script
        if error_data:
            error_type = error_data.get("error_type", "UNKNOWN_ERROR")
            exception_message = error_data.get("exception_message", "")
            traceback_str = error_data.get("traceback", "")

            # Default to a generic code if we can't match a known error_type
            exit_code_label = ERROR_TYPE_TO_EXIT_CODE.get(error_type, "ERROR_SCRIPT_FAILED")

            # Use `.format()` to inject the exception and traceback
            return self.exit_codes[exit_code_label].format(exception=exception_message, traceback=traceback_str)
        # 2) If we reach here, _error.json exists but is empty or doesn't exist at all -> no error recorded
        #    Proceed with parsing results.pickle
        try:
            with self.retrieved.base.repository.open("results.pickle", "rb") as handle:
                results = pickle.load(handle)

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
                            self.logger.warning(
                                f"Found extra results that are not included in the output: {results.keys()}"
                            )

                elif len(top_level_output_list) == 1:
                    # Single top-level output, single result
                    top_level_output_list[0]["value"] = self.serialize_output(results, top_level_output_list[0])
                else:
                    return self.exit_codes.ERROR_RESULT_OUTPUT_MISMATCH

                # Store the outputs
                for output in top_level_output_list:
                    self.out(output["name"], output["value"])

        except OSError:
            return self.exit_codes.ERROR_READING_OUTPUT_FILE
        except ValueError as exception:
            self.logger.error(exception)
            return self.exit_codes.ERROR_INVALID_OUTPUT

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
                        serialized_result[key] = general_serializer(value, serializers=self.serializers)
                return serialized_result
            else:
                self.logger.error(f"Expected a dict for namespace '{name}', got {type(result)}.")
                return self.exit_codes.ERROR_INVALID_OUTPUT
        else:
            return general_serializer(result, serializers=self.serializers)
