import pytest
from src import chemistry
import yaml


@pytest.fixture
def initial_state():
    return 'C'


@pytest.fixture
def atom_types():
    return ["C", "O", "N"]


@pytest.fixture
def allowed_ring_sizes():
    return [3, 4, 5, 6]


@pytest.fixture
def config():
    with open('etc/config.yml', 'r') as file:
        config = yaml.safe_load(file)
    return config


@pytest.fixture
def starting_actions(initial_state, config):
    return chemistry.get_valid_actions(
        initial_state,
        config['chemistry'],
        reset=True
    )


def test_starting_actions(starting_actions):
    assert type(starting_actions) is set
    assert all(x in starting_actions for x in [
               'C#C', 'C=C', 'CC', 'C#N', 'C=N', 'CN', 'C=O', 'CO'])
    assert len(starting_actions) == 8
