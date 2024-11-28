import pytest

pytest_plugins = "aiida.tools.pytest_fixtures"


@pytest.fixture(scope="session", autouse=True)
def aiida_profile(aiida_config, aiida_profile_factory):
    """Create and load a profile with RabbitMQ as broker."""
    with aiida_profile_factory(aiida_config, broker_backend="core.rabbitmq") as profile:
        yield profile


@pytest.fixture
def fixture_localhost(aiida_localhost):
    """Return a localhost `Computer`."""
    localhost = aiida_localhost
    localhost.set_default_mpiprocs_per_machine(1)
    return localhost
