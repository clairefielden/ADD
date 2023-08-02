import rdkit.Chem
import rdkit.Chem.Draw
import rdkit.Chem.QED
from . import chemistry
import torch
import torchvision.transforms as transforms
import os
import errno
from . import docking
import abc
import vina
import meeko

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

transform = transforms.Compose([
    transforms.ToTensor()
])


class RLEnvironmentInterface(abc.ABC):
    def __init__(self, config) -> None:
        self.config = config

    @abc.abstractmethod
    def step(self, action) -> tuple:
        pass

    @abc.abstractmethod
    def reset(self):
        pass


class MoleculeEnvironment(RLEnvironmentInterface):
    def __init__(self, config) -> None:
        super().__init__(config)
        self.state = self.config['training']['starting_molecule']
        self.step_counter = 0
        self.actions_sequence = [self.state]
        self.valid_actions = list(chemistry.get_valid_actions(
            self.state, self.config['chemistry'], reset=True)
        )
        docking.prepare_receptor(self.config['pfpi4k-receptor'])

    @ abc.abstractmethod
    def docking_function():
        pass

    def step(self, action) -> tuple:
        self.step_counter += 1
        terminal = torch.tensor([[1]], device=device)

        if action >= len(self.valid_actions):
            raise ValueError("Invalid Action")

        self.state = self.valid_actions[action]
        self.actions_sequence.append(self.state)
        self.valid_actions = list(chemistry.get_valid_actions(
            self.state, self.config['chemistry']))
        options = torch.stack(
            [
                torch.cat(
                    [
                        torch.tensor(chemistry.get_fingerprint(
                            rdkit.Chem.MolFromSmiles(option_mol)
                        ),
                            dtype=torch.float32,
                            device=device),
                        self.get_steps_remaining()
                    ]
                )
                for option_mol in self.valid_actions
            ]
        )

        reward_tensor = torch.tensor([0], dtype=torch.float32, device=device)

        if ((self.step_counter >= self.config['training']['max_steps']) or len(options) == 0):
            terminal = torch.tensor([[0]], device=device)
            ba = self.docking_function()
            reward_tensor = torch.tensor(
                [-ba], dtype=torch.float32, device=device)

        return options, reward_tensor, terminal

    def reset(self):
        self.state = self.config.training.starting_molecule
        self.step_counter = 0
        self.actions_sequence = [self.state]
        self.valid_actions = list(chemistry.get_valid_actions(
            self.state, self.config['chemistry'], reset=True)
        )

        options = torch.stack(
            [
                torch.cat(
                    [
                        torch.tensor(chemistry.get_fingerprint(
                            rdkit.Chem.MolFromSmiles(option_mol)
                        ),
                            dtype=torch.float32,
                            device=device),
                        self.get_steps_remaining()
                    ]
                )
                for option_mol in self.valid_actions
            ]
        )

        return options

    # other methods
    def get_state_fingerprint(self):
        return torch.tensor(chemistry.get_fingerprint(self.state),
                            dtype=torch.float32,
                            device=device)

    def get_steps_remaining(self):
        steps_remaining = self.config['training']['max_steps'] - \
            self.step_counter
        return torch.tensor([steps_remaining],
                            dtype=torch.float32,
                            device=device)

    def setTarget(self, target):
        self.target = fingerprint_rewards.get_fingerprint(target)

    def setStartingMol(self, mol):
        self.starting_mol = mol
        self.reset()

    def render(self):
        im = rdkit.Chem.Draw.MolToImage(rdkit.Chem.MolFromSmiles(self.state))
        im.show()

    def get_state_image_tensor(self):
        im = rdkit.Chem.Draw.MolToImage(rdkit.Chem.MolFromSmiles(self.state))
        return transform(im)

    def get_options_image_tensor(self):
        im = rdkit.Chem.Draw.MolsToGridImage(
            [rdkit.Chem.MolFromSmiles(smiles)]
            for smiles in self.valid_actions
        )
        return transform(im)

    def render_options(self):
        im = rdkit.Chem.Draw.MolsToGridImage(
            [rdkit.Chem.MolFromSmiles(smiles)]
            for smiles in self.valid_actions
        )
        im.show()

    def save(self, path):
        rdkit.Chem.Draw.MolToFile(rdkit.Chem.MolFromSmiles(self.state), path)

    def saveStateOptions(self, path):
        im = rdkit.Chem.Draw.MolsToGridImage(
            [rdkit.Chem.MolFromSmiles(smiles)]
            for smiles in self.valid_actions
        )
        im.save(path)

    def close(self):
        pass

    def shape_reward(self, reward):
        # shapes rewards by steps remaining
        # steps earlier on in the episode have their reward reduced
        # so that the agent has to obtain high BA states near the end of the
        # episode instead of early in the episode in order to maximize reward
        shaped_reward = reward * (self.config['training']['reward_discount_factor'] **
                                  (self.config['training']['max_steps'] - self.step_counter))
        return shaped_reward

    def qed_reward(self, mol):
        return rdkit.Chem.QED.qed(mol)


class VinaEnvironment(MoleculeEnvironment):
    def __init__(self, config) -> None:
        super().__init__(config)
        self.vina = vina.Vina(sf_name='vina')
        self.vina.set_receptor(self.config['docking-vina']['receptor'])

    def docking_function(self):
        return docking.evaluate_state_vina(
            self.vina,
            self.state,
            self.config['docking-vina']['center'],
            self.config['docking-vina']['box_size'],
            self.config['docking-vina']['exhaustiveness'],
            self.config['docking-vina']['n_poses'])


class AutodockEnvironment(MoleculeEnvironment):
    def __init__(self, config) -> None:
        super().__init__(config)
        self.grid = docking.read_grid_file(
            self.config['docking_autodock_pfpi4k']['grid_file'])
        self.gpu_available = torch.cuda.is_available()

    def docking_function(self):
        return docking.evaluate_state(
            rdkit.Chem.MolFromSmiles(self.state),
            self.config['docking_autodock_pfpi4k']['receptor'],
            self.grid, self.config['docking_autodock_pfpi4k']['nruns'],
            self.gpu_available)
