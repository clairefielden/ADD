import pytest
import rdkit
from vina import Vina
from src import docking
import yaml
import torch


@pytest.fixture
def config():
    with open('etc/config.yml', 'r') as file:
        config = yaml.safe_load(file)
    return config


@pytest.fixture
def vina(config) -> Vina:
    vina = Vina(sf_name='vina')
    vina.set_receptor(config['docking-vina']['receptor'])
    return vina


@pytest.fixture
def target_receptor(config):
    return config['docking_autodock_pfpi4k']['receptor']


@pytest.fixture
def docking_box():
    return 'pfpi4k-docking-box.txt'


@pytest.fixture
def dummy_mol():
    return 'C'


@pytest.fixture
def grid(config):
    return docking.read_grid_file(
        config['docking_autodock_pfpi4k']['grid_file'])


# @pytest.mark.skip
# test initial state autodock docking
def test_evaluate_state(dummy_mol, target_receptor, config, grid):
    assert docking.evaluate_state(
        rdkit.Chem.MolFromSmiles(dummy_mol),
        target_receptor,
        grid,
        config['docking_autodock_pfpi4k']['nruns'],
        torch.cuda.is_available()) == pytest.approx(-2.5, rel=0.1)


# test initial state vina docking
def test_evaluate_state_vina(vina, dummy_mol, config):
    assert docking.evaluate_state_vina(
        vina,
        dummy_mol,
        config['docking-vina']['center'],
        config['docking-vina']['box_size'],
        config['docking-vina']['exhaustiveness'],
        config['docking-vina']['n_poses']
    ) == pytest.approx(-1.2, rel=0.1)
