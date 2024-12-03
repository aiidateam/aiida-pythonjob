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


@patch("aiida_pythonjob.data.serializer.load_config")
@patch("aiida_pythonjob.data.serializer.entry_points")
def test_get_serializer_from_entry_points(mock_entry_points, mock_load_config):
    # Mock the configuration
    mock_load_config.return_value = {
        "serializers": {
            "excludes": ["excluded_entry"],
        }
    }
    # Mock entry points
    mock_ep_1 = create_entry_point("xyz.abc.Abc", "xyz.abc:AbcData", "aiida.data")
    mock_ep_2 = create_entry_point("xyz.abc.Bcd", "xyz.abc:BcdData", "aiida.data")
    mock_ep_3 = create_entry_point("xyz.abc.Cde", "xyz.abc:CdeData", "aiida.data")
    mock_ep_4 = create_entry_point("another_xyz.abc.Cde", "another_xyz.abc:CdeData", "aiida.data")

    mock_entry_points.return_value = create_mock_entry_points([mock_ep_1, mock_ep_2, mock_ep_3, mock_ep_4])

    # Import the function and run
    from aiida_pythonjob.data.serializer import get_serializer_from_entry_points

    with pytest.raises(ValueError, match="Duplicate entry points for abc.Cde"):
        get_serializer_from_entry_points()
    # Mock the configuration
    mock_load_config.return_value = {
        "serializers": {
            "excludes": ["excluded_entry"],
            "abc.Cde": "another_xyz.abc.Cde",
        }
    }
    result = get_serializer_from_entry_points()
    # Assert results
    expected = {
        "abc.Abc": [mock_ep_1],
        "abc.Bcd": [mock_ep_2],
        "abc.Cde": [mock_ep_4],
    }
    assert result == expected
