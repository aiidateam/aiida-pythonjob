from __future__ import annotations

import json
import pathlib
import tempfile
from typing import Any

import cloudpickle as pickle
import pytest
from aiida import orm
from aiida.cmdline.utils.common import get_workchain_report
from aiida.common.links import LinkType
from node_graph.socket_spec import namespace

from aiida_pythonjob.data.deserializer import all_deserializers
from aiida_pythonjob.data.serializer import all_serializers

outputs_spec_with_multiple_sub_specs = namespace(a=Any, b=Any, c=Any)


def create_retrieved_folder(result: dict, error: dict | None = None, output_filename="results.pickle"):
    # Create a retrieved ``FolderData`` node with results
    with tempfile.TemporaryDirectory() as tmpdir:
        dirpath = pathlib.Path(tmpdir)
        with open((dirpath / output_filename), "wb") as handle:
            pickle.dump(result, handle)
        error = error or {}
        with open((dirpath / "_error.json"), "w") as handle:
            json.dump(error, handle)
        folder_data = orm.FolderData(tree=dirpath.absolute())
    return folder_data


def create_process_node(
    result: dict, spec_data: dict, error: dict | None = None, output_filename: str = "results.pickle"
):
    node = orm.CalcJobNode()
    node.set_process_type("aiida.calculations:pythonjob.pythonjob")
    retrieved = create_retrieved_folder(result, error=error, output_filename=output_filename)
    for key, value in spec_data.items():
        node.base.attributes.set(key, value)
    node.base.attributes.set("serializers", all_serializers)
    node.base.attributes.set("deserializers", all_deserializers)
    retrieved.base.links.add_incoming(node, link_type=LinkType.CREATE, link_label="retrieved")
    node.store()
    retrieved.store()
    return node


def create_parser(result, spec_data, error: dict | None = None, output_filename: str = "results.pickle"):
    from aiida_pythonjob.parsers import PythonJobParser

    node = create_process_node(result, spec_data, error=error, output_filename=output_filename)
    parser = PythonJobParser(node=node)
    return parser


def test_tuple_result(fixture_localhost):
    result = (1, 2, 3)
    spec_data = {"outputs_spec": outputs_spec_with_multiple_sub_specs.to_dict()}
    parser = create_parser(result, spec_data)
    exit_code = parser.parse()
    assert exit_code is None
    assert len(parser.outputs) == 3


def test_tuple_result_mismatch(fixture_localhost):
    result = (1, 2)
    spec_data = {"outputs_spec": outputs_spec_with_multiple_sub_specs.to_dict()}
    parser = create_parser(result, spec_data)
    exit_code = parser.parse()
    assert exit_code == parser.exit_codes.ERROR_RESULT_OUTPUT_MISMATCH


def test_dict_result(fixture_localhost):
    result = {"a": 1, "b": 2, "c": 3}
    spec_data = {"outputs_spec": namespace(a=Any, b=Any).to_dict()}
    parser = create_parser(result, spec_data)
    exit_code = parser.parse()
    assert exit_code is None
    assert len(parser.outputs) == 2
    report = get_workchain_report(parser.node, levelname="WARNING")
    assert "Found extra results that are not included in the output: ['c']" in report


def test_dict_result_missing(fixture_localhost):
    result = {"a": 1, "b": 2}
    spec_data = {"outputs_spec": outputs_spec_with_multiple_sub_specs.to_dict()}
    parser = create_parser(result, spec_data)
    exit_code = parser.parse()
    assert exit_code == parser.exit_codes.ERROR_MISSING_OUTPUT


def test_dict_result_as_one_output(fixture_localhost):
    result = {"a": 1, "b": 2, "c": 3}
    spec_data = {"outputs_spec": namespace(result=Any).to_dict()}
    parser = create_parser(result, spec_data)
    exit_code = parser.parse()
    assert exit_code is None
    assert len(parser.outputs) == 1
    assert parser.outputs["result"] == result


def test_dict_result_only_show_one_output(fixture_localhost):
    result = {"a": 1, "b": 2}
    spec_data = {"outputs_spec": namespace(a=Any).to_dict()}
    parser = create_parser(result, spec_data)
    parser.parse()
    assert len(parser.outputs) == 1
    assert parser.outputs["a"] == 1
    report = get_workchain_report(parser.node, levelname="WARNING")
    assert "Found extra results that are not included in the output: ['b']" in report


def test_exit_code(fixture_localhost):
    result = {"a": 1, "exit_code": {"status": 0, "message": ""}}
    spec_data = {"outputs_spec": namespace(a=Any).to_dict()}
    parser = create_parser(result, spec_data)
    exit_code = parser.parse()
    assert exit_code is None
    assert parser.outputs["a"] == 1
    result = {"exit_code": {"status": 1, "message": "error"}}
    spec_data = {"outputs_spec": outputs_spec_with_multiple_sub_specs.to_dict()}
    parser = create_parser(result, spec_data)
    exit_code = parser.parse()
    assert exit_code is not None
    assert exit_code.status == 1
    assert exit_code.message == "error"


def test_no_output_file(fixture_localhost):
    result = {"a": 1, "b": 2, "c": 3}
    spec_data = {"outputs_spec": namespace(result=Any).to_dict()}
    parser = create_parser(result, spec_data, output_filename="not_results.pickle")
    exit_code = parser.parse()
    assert exit_code == parser.exit_codes.ERROR_READING_OUTPUT_FILE


@pytest.mark.parametrize(
    "error_type, status",
    [
        ("IMPORT_CLOUDPICKLE_FAILED", 323),
        ("UNPICKLE_INPUTS_FAILED", 324),
        ("UNPICKLE_FUNCTION_FAILED", 325),
        ("FUNCTION_EXECUTION_FAILED", 326),
        ("PICKLE_RESULTS_FAILED", 327),
    ],
)
def test_run_script_error(error_type, status):
    error = {"error_type": error_type, "exception_message": "error", "traceback": "traceback"}
    result = {"a": 1, "exit_code": {"status": 0, "message": ""}}
    spec_data = {"outputs_spec": namespace(a=Any).to_dict()}
    parser = create_parser(result, spec_data, error=error)
    exit_code = parser.parse()
    assert exit_code is not None
    assert exit_code.status == status
