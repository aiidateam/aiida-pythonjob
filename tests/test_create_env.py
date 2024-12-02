from unittest.mock import MagicMock, patch


def test_create_conda_env():
    computer_name = "test_computer"
    env_name = "test_env"
    pip_packages = ["numpy", "pandas"]
    modules = ["qe"]
    variables = {"TEST_VAR": "test_value"}
    conda_deps = ["scipy"]
    python_version = "3.8"
    shell = "posix"

    # Mock the computer and related objects
    mock_computer = MagicMock()
    mock_computer.label = computer_name
    mock_user = MagicMock()
    mock_user.email = "test_user@test.com"
    mock_authinfo = MagicMock()
    mock_transport = MagicMock()
    mock_scheduler = MagicMock()

    mock_authinfo.get_transport.return_value = mock_transport
    mock_computer.get_authinfo.return_value = mock_authinfo
    mock_authinfo.computer.get_scheduler.return_value = mock_scheduler

    # Mock successful transport behavior
    mock_transport.exec_command_wait.return_value = (
        0,  # retval
        "Environment setup is complete.\n",  # stdout
        "",  # stderr
    )

    # Patch `load_computer` and `User.collection.get_default` to return mocked objects
    with (
        patch("aiida_pythonjob.utils.load_computer", return_value=mock_computer),
        patch("aiida_pythonjob.utils.User.collection.get_default", return_value=mock_user),
    ):
        from aiida_pythonjob.utils import create_conda_env

        success, message = create_conda_env(
            computer=computer_name,
            name=env_name,
            pip=pip_packages,
            conda={"dependencies": conda_deps, "channels": ["conda-forge"]},
            modules=modules,
            variables=variables,
            python_version=python_version,
            shell=shell,
        )

        # Assertions for successful case
        assert success is True
        assert message == "Environment setup is complete."

        # Validate that exec_command_wait was called with the generated script
        mock_transport.exec_command_wait.assert_called_once()
        called_script = mock_transport.exec_command_wait.call_args[0][0]
        assert f"conda create -y -n {env_name} python={python_version}" in called_script
        assert "pip install numpy pandas" in called_script
        assert "conda config --prepend channels" in called_script
        assert "module load qe" in called_script
        assert "export TEST_VAR='test_value'" in called_script


def test_create_conda_env_error_handling():
    computer_name = "test_computer"
    env_name = "test_env"

    # Mock the computer and related objects
    mock_computer = MagicMock()
    mock_computer.label = computer_name
    mock_user = MagicMock()
    mock_user.email = "test_user@test.com"
    mock_authinfo = MagicMock()
    mock_transport = MagicMock()
    mock_scheduler = MagicMock()

    # Mock error in transport
    mock_transport.exec_command_wait.return_value = (
        1,  # retval
        "",  # stdout
        "Error creating environment",  # stderr
    )

    mock_authinfo.get_transport.return_value = mock_transport
    mock_authinfo.computer.get_scheduler.return_value = mock_scheduler
    mock_computer.get_authinfo.return_value = mock_authinfo

    # Patch `load_computer` and `User.collection.get_default` to return mocked objects
    with (
        patch("aiida_pythonjob.utils.load_computer", return_value=mock_computer),
        patch("aiida_pythonjob.utils.User.collection.get_default", return_value=mock_user),
    ):
        from aiida_pythonjob.utils import create_conda_env

        success, message = create_conda_env(
            computer=computer_name,
            name=env_name,
        )

        # Assertions for failure case
        assert success is False
        assert "The command returned a non-zero return code" in message
