import sys
from importlib.metadata import EntryPoint
from unittest.mock import patch

import pytest


# Helper function to mock EntryPoint creation
def create_entry_point(name, value, group):
    return EntryPoint(name=name, value=value, group=group)


def create_mock_entry_points(entry_point_list):
    if sys.version_info >= (3, 10):
        # Mock the EntryPoints object for Python 3.10+
        # Conditional import for EntryPoints
        from importlib.metadata import EntryPoints

        return EntryPoints(entry_point_list)
    else:
        # Return a dictionary for older Python versions
        return {"aiida.data": entry_point_list}


@patch("aiida_pythonjob.data.serializer.entry_points")
def test_get_serializers(mock_entry_points):
    # Mock the configuration
    from aiida_pythonjob.config import config

    config["serializers"] = {}
    # Mock entry points
    mock_ep_1 = create_entry_point("xyz.abc.Abc", "xyz.abc:AbcData", "aiida.data")
    mock_ep_2 = create_entry_point("xyz.abc.Bcd", "xyz.abc:BcdData", "aiida.data")
    mock_ep_3 = create_entry_point("xyz.abc.Cde", "xyz.abc:CdeData", "aiida.data")
    mock_ep_4 = create_entry_point("another_xyz.abc.Cde", "another_xyz.abc:CdeData", "aiida.data")

    mock_entry_points.return_value = create_mock_entry_points([mock_ep_1, mock_ep_2, mock_ep_3, mock_ep_4])

    # Import the function and run
    from aiida_pythonjob.data.serializer import get_serializers

    with pytest.raises(ValueError, match="Duplicate entry points for abc.Cde"):
        get_serializers()
    # Mock the configuration
    config["serializers"] = {
        "abc.Cde": "another_xyz.abc.CdeData",
    }
    result = get_serializers()
    # Assert results
    expected = {"abc.Abc": "xyz.abc.AbcData", "abc.Bcd": "xyz.abc.BcdData", "abc.Cde": "another_xyz.abc.CdeData"}
    print("result", result)
    assert result == expected
