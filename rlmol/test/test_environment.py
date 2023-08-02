import pytest
from src import environment, docking
import torch
import random
import yaml


@pytest.fixture
def config():
    with open('etc/config.yml', 'r') as file:
        config = yaml.safe_load(file)
    return config


@pytest.fixture
def target_receptor():
    return 'pfpi4k-no-ligand'


@pytest.fixture
def docking_box():
    return 'pfpi4k-docking-box.txt'


@pytest.fixture
def basic_environment(config):
    return environment.AutodockEnvironment(config)


@pytest.fixture
def vina_environment(config):
    return environment.VinaEnvironment(config)


# @pytest.mark.skip
def test_initialise(basic_environment):
    assert basic_environment.state == 'C'
    assert basic_environment.actions_sequence == ['C']
    assert all(x in basic_environment.valid_actions for x in [
               'C#C', 'C=C', 'CC', 'C#N', 'C=N', 'CN', 'C=O', 'CO'])
    assert len(basic_environment.valid_actions) == 8
    # CH4 should be approx -2.5
    assert basic_environment.docking_function() == pytest.approx(-2.5, rel=0.1)


# @pytest.mark.skip
def test_initialise_vina_environemnt(vina_environment):
    assert vina_environment.state == 'C'
    assert vina_environment.actions_sequence == ['C']
    assert all(x in vina_environment.valid_actions for x in [
               'C#C', 'C=C', 'CC', 'C#N', 'C=N', 'CN', 'C=O', 'CO'])
    assert len(vina_environment.valid_actions) == 8
    # CH4 should be approx -2.5
    assert vina_environment.docking_function() == pytest.approx(-1.2, rel=0.1)


# @pytest.mark.skip
def test_step_vina(vina_environment):
    chosen_state = vina_environment.valid_actions[2]
    options, reward, terminal = vina_environment.step(2)
    post_step_state = vina_environment.state
    assert post_step_state == chosen_state  # should move into chosen state
    assert type(options) == torch.Tensor
    assert type(reward) == torch.Tensor
    assert type(terminal) == torch.Tensor
    assert reward.size() == torch.Size([1])
    assert reward.dim() == 1
    assert terminal.size() == torch.Size([1, 1])  # contains 1 row
    assert terminal.dim() == 2  # dim is 2 because flag is store in its own row
    assert reward.item() == 0  # reward should be 0 for non-terminal states
    assert terminal.item() == 1  # terminal flag should be 1 for non-terminal states
    assert options.dim() == 2  # number of tensor dimensions
    assert options.size(dim=1) == 2049  # 2048 bits + 1 int for steps remaining
    assert options[0][2048].item() == 39  # steps remaining after 1 action


# @pytest.mark.skip
def test_step(basic_environment):
    chosen_state = basic_environment.valid_actions[2]
    options, reward, terminal = basic_environment.step(2)
    post_step_state = basic_environment.state
    assert post_step_state == chosen_state  # should move into chosen state
    assert type(options) == torch.Tensor
    assert type(reward) == torch.Tensor
    assert type(terminal) == torch.Tensor
    assert reward.size() == torch.Size([1])
    assert reward.dim() == 1
    assert terminal.size() == torch.Size([1, 1])  # contains 1 row
    assert terminal.dim() == 2  # dim is 2 because flag is store in its own row
    assert reward.item() == 0  # reward should be 0 for non-terminal states
    assert terminal.item() == 1  # terminal flag should be 1 for non-terminal states
    assert options.dim() == 2  # number of tensor dimensions
    assert options.size(dim=1) == 2049  # 2048 bits + 1 int for steps remaining
    assert options[0][2048].item() == 39  # steps remaining after 1 action


# @pytest.mark.skip
def test_episode_vina(vina_environment):
    # intial step
    action = random.randrange(len(vina_environment.valid_actions))
    assert vina_environment.step_counter == 0
    # intervening steps
    for i in range(1, vina_environment.config['training']['max_steps']):
        options, reward, terminal = vina_environment.step(action)
        assert vina_environment.step_counter == i
        action = random.randrange(options.size(dim=0))
        state = options[action]
        steps_remaining = state[2048]
        assert steps_remaining == vina_environment.config['training']['max_steps'] - \
            vina_environment.step_counter
        assert terminal.item() == 1
    # before final step
    assert state[2048] == 1  # 1 step remaining
    assert terminal.item() == 1  # terminal flag not set
    # final step - should return reward, steps remaining should be 0
    options, reward, terminal = vina_environment.step(action)
    action = random.randrange(len(vina_environment.valid_actions))
    state = options[action]
    assert state[2048] == 0  # no steps remaining
    assert terminal.item() == 0  # terminal flag set
    # assert reward == pytest.approx(2.5, rel=0.1)


# @pytest.mark.skip
def test_episode(basic_environment):
    # intial step
    action = random.randrange(len(basic_environment.valid_actions))
    assert basic_environment.step_counter == 0
    # intervening steps
    for i in range(1, basic_environment.config['training']['max_steps']):
        options, reward, terminal = basic_environment.step(action)
        assert basic_environment.step_counter == i
        action = random.randrange(options.size(dim=0))
        state = options[action]
        steps_remaining = state[2048]
        assert steps_remaining == basic_environment.config['training']['max_steps'] - \
            basic_environment.step_counter
        assert terminal.item() == 1
    # before final step
    assert state[2048] == 1  # 1 step remaining
    assert terminal.item() == 1  # terminal flag not set
    # final step - should return reward, steps remaining should be 0
    options, reward, terminal = basic_environment.step(action)
    action = random.randrange(len(basic_environment.valid_actions))
    state = options[action]
    assert state[2048] == 0  # no steps remaining
    assert terminal.item() == 0  # terminal flag set
    # assert reward == pytest.approx(2.5, rel=0.1)
