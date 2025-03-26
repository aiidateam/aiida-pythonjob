"""Calcjob to run a Python function on a remote computer, either via raw source code or a pickled function."""

from __future__ import annotations

import pathlib
import typing as t

from aiida.common.datastructures import CalcInfo, CodeInfo
from aiida.common.folders import Folder
from aiida.common.lang import override
from aiida.engine import CalcJob, CalcJobProcessSpec
from aiida.orm import (
    Data,
    Dict,
    FolderData,
    List,
    RemoteData,
    SinglefileData,
    Str,
    to_aiida_type,
)

__all__ = ("PythonJob",)


class PythonJob(CalcJob):
    """Calcjob to run a Python function on a remote computer.

    Supports two modes:
    1) Loading a pickled function object (function_data.pickled_function).
    2) Embedding raw source code for the function (function_data.source_code).
    """

    _internal_retrieve_list = []
    _retrieve_singlefile_list = []
    _retrieve_temporary_list = []

    _DEFAULT_INPUT_FILE = "script.py"
    _DEFAULT_OUTPUT_FILE = "aiida.out"
    _DEFAULT_PARENT_FOLDER_NAME = "./parent_folder/"
    _SOURCE_CODE_KEY = "source_code"

    @classmethod
    def define(cls, spec: CalcJobProcessSpec) -> None:  # type: ignore[override]
        """Define the process specification, including its inputs, outputs and known exit codes."""
        super().define(spec)
        spec.input_namespace("function_data", dynamic=True, required=True)
        spec.input("function_data.input_ports", valid_type=Dict, serializer=to_aiida_type, required=False)
        spec.input("function_data.output_ports", valid_type=Dict, serializer=to_aiida_type, required=False)
        spec.input("process_label", valid_type=Str, serializer=to_aiida_type, required=False)
        spec.input_namespace("function_inputs", valid_type=Data, required=False)
        spec.input(
            "parent_folder",
            valid_type=(RemoteData, FolderData, SinglefileData),
            required=False,
            help="Use a local or remote folder as parent folder (for restarts and similar)",
        )
        spec.input(
            "parent_folder_name",
            valid_type=Str,
            required=False,
            serializer=to_aiida_type,
            help="""Default name of the subfolder to create in the working directory
            where the files from parent_folder are placed.""",
        )
        spec.input(
            "parent_output_folder",
            valid_type=Str,
            default=None,
            required=False,
            serializer=to_aiida_type,
            help="Name of the subfolder inside 'parent_folder' from which you want to copy the files",
        )
        spec.input_namespace(
            "upload_files",
            valid_type=(FolderData, SinglefileData),
            required=False,
            help="The folder/files to upload",
        )
        spec.input_namespace(
            "copy_files",
            valid_type=(RemoteData,),
            required=False,
            help="The folder/files to copy from the remote computer",
        )
        spec.input(
            "additional_retrieve_list",
            valid_type=List,
            default=None,
            required=False,
            serializer=to_aiida_type,
            help="Additional filenames to retrieve from the remote work directory",
        )
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
        spec.outputs.dynamic = True
        # set default options (optional)
        spec.inputs["metadata"]["options"]["parser_name"].default = "pythonjob.pythonjob"
        spec.inputs["metadata"]["options"]["input_filename"].default = "script.py"
        spec.inputs["metadata"]["options"]["output_filename"].default = "aiida.out"
        spec.inputs["metadata"]["options"]["resources"].default = {
            "num_machines": 1,
            "num_mpiprocs_per_machine": 1,
        }
        spec.exit_code(
            310,
            "ERROR_READING_OUTPUT_FILE",
            invalidates_cache=True,
            message="The output file could not be read.",
        )
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
            322,
            "ERROR_IMPORT_CLOUDPICKLE_FAILED",
            invalidates_cache=True,
            message="Importing cloudpickle failed.\n{exception}\n{traceback}",
        )
        spec.exit_code(
            323,
            "ERROR_UNPICKLE_INPUTS_FAILED",
            invalidates_cache=True,
            message="Failed to unpickle inputs.\n{exception}\n{traceback}",
        )
        spec.exit_code(
            324,
            "ERROR_UNPICKLE_FUNCTION_FAILED",
            invalidates_cache=True,
            message="Failed to unpickle user function.\n{exception}\n{traceback}",
        )
        spec.exit_code(
            325,
            "ERROR_FUNCTION_EXECUTION_FAILED",
            invalidates_cache=True,
            message="Function execution failed.\n{exception}\n{traceback}",
        )
        spec.exit_code(
            326,
            "ERROR_PICKLE_RESULTS_FAILED",
            invalidates_cache=True,
            message="Failed to pickle results.\n{exception}\n{traceback}",
        )
        spec.exit_code(
            327,
            "ERROR_SCRIPT_FAILED",
            invalidates_cache=True,
            message="The script failed for an unknown reason.\n{exception}\n{traceback}",
        )

    def get_function_name(self) -> str:
        """Return the name of the function to run."""
        if "name" in self.inputs.function_data:
            name = self.inputs.function_data.name
        else:
            name = "anonymous_function"
        return name

    def _build_process_label(self) -> str:
        """Use the function name or an explicit label as the process label."""
        if "process_label" in self.inputs:
            return self.inputs.process_label.value
        else:
            name = self.get_function_name()
            return f"PythonJob<{name}>"

    @override
    def _setup_db_record(self) -> None:
        """Set up the database record for the process."""
        super()._setup_db_record()
        if "source_code" in self.inputs.function_data:
            self.node.base.attributes.set(self._SOURCE_CODE_KEY, self.inputs.function_data.source_code)

    def on_create(self) -> None:
        """Called when a Process is created."""
        super().on_create()
        self.node.label = self._build_process_label()

    def prepare_for_submission(self, folder: Folder) -> CalcInfo:
        """Prepare the calculation for submission.

        1) Write the python script to the folder, depending on the mode (source vs. pickled function).
        2) Write the inputs to a pickle file and save it to the folder.
        """
        import cloudpickle as pickle

        from aiida_pythonjob.calculations.utils import generate_script_py
        from aiida_pythonjob.data.deserializer import deserialize_to_raw_python_data

        dirpath = pathlib.Path(folder._abspath)

        # Prepare the dictionary of input arguments for the function
        inputs: dict[str, t.Any]
        if self.inputs.function_inputs:
            inputs = dict(self.inputs.function_inputs)
        else:
            inputs = {}

        # Prepare the final subfolder name for the parent folder
        if "parent_folder_name" in self.inputs:
            parent_folder_name = self.inputs.parent_folder_name.value
        else:
            parent_folder_name = self._DEFAULT_PARENT_FOLDER_NAME

        # Build the Python script
        source_code = self.node.base.attributes.get(self._SOURCE_CODE_KEY, None)
        pickled_function = self.inputs.function_data.pickled_function
        function_name = self.get_function_name()  # or some user-defined name
        script_content = generate_script_py(
            pickled_function=pickled_function,
            source_code=source_code,
            function_name=function_name,
        )

        # Write the script to the working folder
        with folder.open(self.options.input_filename, "w", encoding="utf8") as handle:
            handle.write(script_content)

        # Symlink or copy approach for the parent folder
        symlink = True
        remote_copy_list = []
        local_copy_list = []
        remote_symlink_list = []
        remote_list = remote_symlink_list if symlink else remote_copy_list

        source = self.inputs.get("parent_folder", None)
        if source is not None:
            if isinstance(source, RemoteData):
                # Possibly append parent_output_folder path
                dirpath_remote = pathlib.Path(source.get_remote_path())
                if self.inputs.parent_output_folder is not None:
                    dirpath_remote /= self.inputs.parent_output_folder.value
                remote_list.append(
                    (
                        source.computer.uuid,
                        str(dirpath_remote),
                        parent_folder_name,
                    )
                )
            elif isinstance(source, FolderData):
                dirname = self.inputs.parent_output_folder.value if self.inputs.parent_output_folder is not None else ""
                local_copy_list.append((source.uuid, dirname, parent_folder_name))
            elif isinstance(source, SinglefileData):
                local_copy_list.append((source.uuid, source.filename, source.filename))

        # Upload additional files
        if "upload_files" in self.inputs:
            upload_files = self.inputs.upload_files
            for key, src in upload_files.items():
                # replace "_dot_" with "." in the key
                new_key = key.replace("_dot_", ".")
                if isinstance(src, FolderData):
                    local_copy_list.append((src.uuid, "", new_key))
                elif isinstance(src, SinglefileData):
                    local_copy_list.append((src.uuid, src.filename, src.filename))
                else:
                    raise ValueError(
                        f"Input file/folder '{key}' of type {type(src)} is not supported. "
                        "Only AiiDA SinglefileData and FolderData are allowed."
                    )

        # Copy remote data if any
        if "copy_files" in self.inputs:
            copy_files = self.inputs.copy_files
            for key, src in copy_files.items():
                new_key = key.replace("_dot_", ".")
                dirpath_remote = pathlib.Path(src.get_remote_path())
                remote_list.append((src.computer.uuid, str(dirpath_remote), new_key))

        # Create a pickle file for the user input values
        input_values = {}
        if "deserializers" in self.inputs and self.inputs.deserializers:
            deserializers = self.inputs.deserializers.get_dict()
            # replace "__dot__" with "." in the keys
            deserializers = {k.replace("__dot__", "."): v for k, v in deserializers.items()}
        else:
            deserializers = None
        input_values = deserialize_to_raw_python_data(inputs, deserializers=deserializers)

        filename = "inputs.pickle"
        with folder.open(filename, "wb") as handle:
            pickle.dump(input_values, handle)

        # If using a pickled function, we also need to upload the function pickle
        if pickled_function:
            # create a SinglefileData object for the pickled function
            function_pkl_fname = "function.pkl"
            with folder.open(function_pkl_fname, "wb") as handle:
                handle.write(pickled_function)

        # create a singlefiledata object for the pickled data
        file_data = SinglefileData(file=f"{dirpath}/{filename}")
        file_data.store()
        local_copy_list.append((file_data.uuid, file_data.filename, filename))

        codeinfo = CodeInfo()
        codeinfo.stdin_name = self.options.input_filename
        codeinfo.stdout_name = self.options.output_filename
        codeinfo.code_uuid = self.inputs.code.uuid

        calcinfo = CalcInfo()
        calcinfo.codes_info = [codeinfo]
        calcinfo.local_copy_list = local_copy_list
        calcinfo.remote_copy_list = remote_copy_list
        calcinfo.remote_symlink_list = remote_symlink_list
        calcinfo.retrieve_list = ["results.pickle", self.options.output_filename, "_error.json"]
        if self.inputs.additional_retrieve_list is not None:
            calcinfo.retrieve_list += self.inputs.additional_retrieve_list.get_list()
        calcinfo.retrieve_list += self._internal_retrieve_list

        calcinfo.retrieve_temporary_list = self._retrieve_temporary_list
        calcinfo.retrieve_singlefile_list = self._retrieve_singlefile_list

        return calcinfo
