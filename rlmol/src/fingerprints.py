from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import QED


def get_fingerprint(molecule):
    """Gets the morgan fingerprint of the target molecule.

    Args:
        molecule: Chem.Mol. The current molecule.

    Returns:
        rdkit.ExplicitBitVect. The fingerprint of the target.
    """
    return AllChem.GetMorganFingerprint(molecule, radius=2)


def get_similarity(current_mol, target_mol_fingerprint):
    """Gets the similarity between the current molecule and the target molecule.

    Args:
        smiles: String. The SMILES string for the current molecule.

    Returns:
        Float. The Tanimoto similarity.
    """
    current_mol_fingerprint = get_fingerprint(current_mol)

    return DataStructs.TanimotoSimilarity(target_mol_fingerprint,
                                          current_mol_fingerprint)


def rewardv0(current_similarity):
    if current_similarity > 0:
        return current_similarity
    else:
        return -1


def reward(current_mol_smiles, target_mol_smiles, max_steps, counter,
           similarity_weight=0.5, discount_factor=0.999):
    """Calculates the reward of the current state.

    The reward is defined as a tuple of the similarity and QED value.

    Returns:
        A tuple of the similarity and qed value
    """
    # calculate similarity.
    # if the current molecule does not contain the scaffold of the target,
    # similarity is zero.
    if current_mol_smiles is None:
        return 0.0
    mol = Chem.MolFromSmiles(current_mol_smiles)
    if mol is None:
        return 0.0
    similarity_score = get_similarity(current_mol_smiles, target_mol_smiles)
    # calculate QED
    qed_value = QED.qed(mol)
    reward = (
        similarity_score * similarity_weight +
        qed_value * (1 - similarity_weight))
    discount = discount_factor**(max_steps - counter)
    return reward * discount
