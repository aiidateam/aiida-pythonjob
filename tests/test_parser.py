from __future__ import annotations

import json
import pathlib
import tempfile

import cloudpickle as pickle
import pytest
from aiida import orm
from aiida.cmdline.utils.common import get_workchain_report
from aiida.common.links import LinkType

output_ports_with_multiple_sub_ports = {
    "name": "outputs",
    "identifier": "namespace",
    "ports": [
        {"name": "a", "identifier": "any"},
        {"name": "b", "identifier": "any"},
        {"name": "c", "identifier": "any"},
    ],
}


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
    result: dict, function_data: dict, error: dict | None = None, output_filename: str = "results.pickle"
):
    node = orm.CalcJobNode()
    node.set_process_type("aiida.calculations:pythonjob.pythonjob")
    retrieved = create_retrieved_folder(result, error=error, output_filename=output_filename)
    for key, value in function_data.items():
        node.base.links.add_incoming(value, link_type=LinkType.INPUT_CALC, link_label=f"function_data__{key}")
        value.store()
    retrieved.base.links.add_incoming(node, link_type=LinkType.CREATE, link_label="retrieved")
    node.store()
    retrieved.store()
    return node


def create_parser(result, function_data, error: dict | None = None, output_filename: str = "results.pickle"):
    from aiida_pythonjob.parsers import PythonJobParser

    node = create_process_node(result, function_data, error=error, output_filename=output_filename)
    parser = PythonJobParser(node=node)
    return parser


def test_tuple_result(fixture_localhost):
    result = (1, 2, 3)
    function_data = {"output_ports": orm.Dict(output_ports_with_multiple_sub_ports)}
    parser = create_parser(result, function_data)
    exit_code = parser.parse()
    assert exit_code is None
    assert len(parser.outputs) == 3


def test_tuple_result_mismatch(fixture_localhost):
    result = (1, 2)
    function_data = {"output_ports": orm.Dict(output_ports_with_multiple_sub_ports)}
    parser = create_parser(result, function_data)
    exit_code = parser.parse()
    assert exit_code == parser.exit_codes.ERROR_RESULT_OUTPUT_MISMATCH


def test_dict_result(fixture_localhost):
    result = {"a": 1, "b": 2, "c": 3}
    function_data = {
        "output_ports": orm.Dict(
            {
                "name": "outputs",
                "identifier": "namespace",
                "ports": [
                    {"name": "a", "identifier": "any"},
                    {"name": "b", "identifier": "any"},
                ],
            }
        )
    }
    parser = create_parser(result, function_data)
    exit_code = parser.parse()
    assert exit_code is None
    assert len(parser.outputs) == 2
    report = get_workchain_report(parser.node, levelname="WARNING")
    assert "Found extra results that are not included in the output: dict_keys(['c'])" in report


def test_dict_result_missing(fixture_localhost):
    result = {"a": 1, "b": 2}
    function_data = {"output_ports": orm.Dict(output_ports_with_multiple_sub_ports)}
    parser = create_parser(result, function_data)
    exit_code = parser.parse()
    assert exit_code == parser.exit_codes.ERROR_MISSING_OUTPUT


def test_dict_result_as_one_output(fixture_localhost):
    result = {"a": 1, "b": 2, "c": 3}
    function_data = {
        "output_ports": orm.Dict(
            {"name": "outputs", "identifier": "namespace", "ports": [{"name": "result", "identifier": "any"}]}
        )
    }
    parser = create_parser(result, function_data)
    exit_code = parser.parse()
    assert exit_code is None
    assert len(parser.outputs) == 1
    assert parser.outputs["result"] == result


def test_dict_result_only_show_one_output(fixture_localhost):
    result = {"a": 1, "b": 2}
    function_data = {
        "output_ports": orm.Dict(
            {"name": "outputs", "identifier": "namespace", "ports": [{"name": "a", "identifier": "any"}]}
        )
    }
    parser = create_parser(result, function_data)
    parser.parse()
    assert len(parser.outputs) == 1
    assert parser.outputs["a"] == 1
    report = get_workchain_report(parser.node, levelname="WARNING")
    assert "Found extra results that are not included in the output: dict_keys(['b'])" in report


def test_exit_code(fixture_localhost):
    result = {"a": 1, "exit_code": {"status": 0, "message": ""}}
    function_data = {
        "output_ports": orm.Dict(
            {"name": "outputs", "identifier": "namespace", "ports": [{"name": "a", "identifier": "any"}]}
        )
    }
    parser = create_parser(result, function_data)
    exit_code = parser.parse()
    assert exit_code is None
    assert parser.outputs["a"] == 1
    #
    result = {"exit_code": {"status": 1, "message": "error"}}
    function_data = {"output_ports": orm.Dict(output_ports_with_multiple_sub_ports)}
    parser = create_parser(result, function_data)
    exit_code = parser.parse()
    assert exit_code is not None
    assert exit_code.status == 1
    assert exit_code.message == "error"


def test_no_output_file(fixture_localhost):
    result = {"a": 1, "b": 2, "c": 3}
    function_data = {
        "output_ports": orm.Dict(
            {"name": "outputs", "identifier": "namespace", "ports": [{"name": "result", "identifier": "any"}]}
        )
    }
    parser = create_parser(result, function_data, output_filename="not_results.pickle")
    exit_code = parser.parse()
    assert exit_code == parser.exit_codes.ERROR_READING_OUTPUT_FILE


@pytest.mark.parametrize(
    "error_type, status",
    [
        ("IMPORT_CLOUDPICKLE_FAILED", 322),
        ("UNPICKLE_INPUTS_FAILED", 323),
        ("UNPICKLE_FUNCTION_FAILED", 324),
        ("FUNCTION_EXECUTION_FAILED", 325),
        ("PICKLE_RESULTS_FAILED", 326),
    ],
)
def test_run_script_error(error_type, status):
    error = {"error_type": error_type, "exception_message": "error", "traceback": "traceback"}
    result = {"a": 1, "exit_code": {"status": 0, "message": ""}}
    function_data = {
        "output_ports": orm.Dict(
            {"name": "outputs", "identifier": "namespace", "ports": [{"name": "a", "identifier": "any"}]}
        )
    }
    parser = create_parser(result, function_data, error=error)
    exit_code = parser.parse()
    assert exit_code is not None
    assert exit_code.status == status
