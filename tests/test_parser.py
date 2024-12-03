import pathlib
import tempfile

import cloudpickle as pickle
from aiida import orm
from aiida.cmdline.utils.common import get_workchain_report
from aiida.common.links import LinkType
from aiida_pythonjob.parsers import PythonJobParser


def create_retrieved_folder(result: dict, output_filename="results.pickle"):
    # Create a retrieved ``FolderData`` node with results
    with tempfile.TemporaryDirectory() as tmpdir:
        dirpath = pathlib.Path(tmpdir)
        with open((dirpath / output_filename), "wb") as handle:
            pickle.dump(result, handle)
        folder_data = orm.FolderData(tree=dirpath.absolute())
    return folder_data


def create_process_node(result: dict, function_data: dict, output_filename: str = "results.pickle"):
    node = orm.CalcJobNode()
    node.set_process_type("aiida.calculations:pythonjob.pythonjob")
    function_data = orm.Dict(function_data)
    retrieved = create_retrieved_folder(result, output_filename=output_filename)
    node.base.links.add_incoming(function_data, link_type=LinkType.INPUT_CALC, link_label="function_data")
    retrieved.base.links.add_incoming(node, link_type=LinkType.CREATE, link_label="retrieved")
    function_data.store()
    node.store()
    retrieved.store()
    return node


def create_parser(result, function_data, output_filename="results.pickle"):
    node = create_process_node(result, function_data, output_filename=output_filename)
    parser = PythonJobParser(node=node)
    return parser


def test_tuple_result(fixture_localhost):
    result = (1, 2, 3)
    function_data = {"outputs": [{"name": "a"}, {"name": "b"}, {"name": "c"}]}
    parser = create_parser(result, function_data)
    exit_code = parser.parse()
    assert exit_code is None
    assert len(parser.outputs) == 3


def test_tuple_result_mismatch(fixture_localhost):
    result = (1, 2)
    function_data = {"outputs": [{"name": "a"}, {"name": "b"}, {"name": "c"}]}
    parser = create_parser(result, function_data)
    exit_code = parser.parse()
    assert exit_code == parser.exit_codes.ERROR_RESULT_OUTPUT_MISMATCH


def test_dict_result(fixture_localhost):
    result = {"a": 1, "b": 2, "c": 3}
    function_data = {"outputs": [{"name": "a"}, {"name": "b"}]}
    parser = create_parser(result, function_data)
    exit_code = parser.parse()
    assert exit_code is None
    assert len(parser.outputs) == 2
    report = get_workchain_report(parser.node, levelname="WARNING")
    assert "Found extra results that are not included in the output: dict_keys(['c'])" in report


def test_dict_result_missing(fixture_localhost):
    result = {"a": 1, "b": 2}
    function_data = {"outputs": [{"name": "a"}, {"name": "b"}, {"name": "c"}]}
    parser = create_parser(result, function_data)
    exit_code = parser.parse()
    assert exit_code == parser.exit_codes.ERROR_MISSING_OUTPUT


def test_dict_result_as_one_output(fixture_localhost):
    result = {"a": 1, "b": 2, "c": 3}
    function_data = {"outputs": [{"name": "result"}]}
    parser = create_parser(result, function_data)
    exit_code = parser.parse()
    assert exit_code is None
    assert len(parser.outputs) == 1
    assert parser.outputs["result"] == result


def test_dict_result_only_show_one_output(fixture_localhost):
    result = {"a": 1, "b": 2}
    function_data = {"outputs": [{"name": "a"}]}
    parser = create_parser(result, function_data)
    parser.parse()
    assert len(parser.outputs) == 1
    assert parser.outputs["a"] == 1
    report = get_workchain_report(parser.node, levelname="WARNING")
    assert "Found extra results that are not included in the output: dict_keys(['b'])" in report


def test_exit_code(fixture_localhost):
    result = {"exit_code": {"status": 1, "message": "error"}}
    function_data = {"outputs": [{"name": "a"}, {"name": "b"}, {"name": "c"}]}
    parser = create_parser(result, function_data)
    exit_code = parser.parse()
    assert exit_code is not None
    assert exit_code.status == 1
    assert exit_code.message == "error"


def test_no_output_file(fixture_localhost):
    result = {"a": 1, "b": 2, "c": 3}
    function_data = {"outputs": [{"name": "result"}]}
    parser = create_parser(result, function_data, output_filename="not_results.pickle")
    exit_code = parser.parse()
    assert exit_code == parser.exit_codes.ERROR_READING_OUTPUT_FILE
