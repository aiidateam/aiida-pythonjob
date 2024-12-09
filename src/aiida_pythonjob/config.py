import json

from aiida.manage.configuration.settings import AiiDAConfigDir


def load_config() -> dict:
    """Load the configuration from the config file."""
    config_file_path = AiiDAConfigDir.get() / "pythonjob.json"
    try:
        with config_file_path.open("r") as f:
            config = json.load(f)
    except FileNotFoundError:
        config = {}
    return config
