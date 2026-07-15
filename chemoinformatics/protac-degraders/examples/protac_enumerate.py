#!/usr/bin/env python3
'''Enumerate connected PROTAC-like molecules from explicit attachment points.'''
# Reference: RDKit 2024.09+ | Verify API if version differs

from rdkit import Chem


# Linker dummy map 1 joins the target fragment; map 2 joins the E3 fragment.
LINKERS = {
    'short_alkyl': '[*:1]CC[*:2]',
    'short_pegylated': '[*:1]COC[*:2]',
    'piperazine_short': '[*:1]N1CCN([*:2])CC1',
    'medium_alkyl': '[*:1]CCCCC[*:2]',
    'medium_pegylated': '[*:1]CCOCC[*:2]',
    'piperazine_acid': '[*:1]N1CCN(CC1)C(=O)[*:2]',
    'triazole': '[*:1]Cn1cc(C[*:2])nn1',
    'long_alkyl': '[*:1]CCCCCCCCCC[*:2]',
    'long_pegylated': '[*:1]COCCOCCOC[*:2]',
    'rigid_piperidine': '[*:1]C(=O)CN1CCC(CCC[*:2])CC1',
    'triazole_extended': '[*:1]CCCCn1cc(CCCC[*:2])nn1',
}


def _dummy_indices(mol, atom_map=None):
    indices = [
        atom.GetIdx() for atom in mol.GetAtoms()
        if atom.GetAtomicNum() == 0
        and (atom_map is None or atom.GetAtomMapNum() == atom_map)
    ]
    return indices


def _connect_dummies(left, right, left_dummy_idx, right_dummy_idx):
    '''Remove two degree-one dummy atoms and bond their neighboring atoms.'''
    left_dummy = left.GetAtomWithIdx(left_dummy_idx)
    right_dummy = right.GetAtomWithIdx(right_dummy_idx)
    if left_dummy.GetDegree() != 1 or right_dummy.GetDegree() != 1:
        raise ValueError('Each attachment dummy must have exactly one neighbor')

    left_neighbor = left_dummy.GetNeighbors()[0].GetIdx()
    right_neighbor = right_dummy.GetNeighbors()[0].GetIdx()
    left_bond = left.GetBondBetweenAtoms(left_dummy_idx, left_neighbor)
    right_bond = right.GetBondBetweenAtoms(right_dummy_idx, right_neighbor)
    if (
        left_bond.GetBondType() != Chem.BondType.SINGLE
        or right_bond.GetBondType() != Chem.BondType.SINGLE
    ):
        raise ValueError('Attachment dummy bonds must be single bonds')
    offset = left.GetNumAtoms()
    combined = Chem.RWMol(Chem.CombineMols(left, right))
    combined.AddBond(left_neighbor, offset + right_neighbor, Chem.BondType.SINGLE)

    for idx in sorted((left_dummy_idx, offset + right_dummy_idx), reverse=True):
        combined.RemoveAtom(idx)
    product = combined.GetMol()
    Chem.SanitizeMol(product)
    return product


def build_protac(target_fragment_smi, linker_smi, e3_fragment_smi):
    '''Connect target-[*:1]linker[*:2]-E3 using explicit dummy atoms.'''
    target = Chem.MolFromSmiles(target_fragment_smi)
    linker = Chem.MolFromSmiles(linker_smi)
    e3 = Chem.MolFromSmiles(e3_fragment_smi)
    if target is None or linker is None or e3 is None:
        raise ValueError('Target, linker, and E3 inputs must be valid SMILES')

    target_dummies = _dummy_indices(target)
    e3_dummies = _dummy_indices(e3)
    linker_target = _dummy_indices(linker, atom_map=1)
    linker_e3 = _dummy_indices(linker, atom_map=2)
    if len(target_dummies) != 1 or len(e3_dummies) != 1:
        raise ValueError('Target and E3 fragments must each contain exactly one [*]')
    if len(linker_target) != 1 or len(linker_e3) != 1:
        raise ValueError('Linker must contain one [*:1] and one [*:2]')

    target_linker = _connect_dummies(
        target, linker, target_dummies[0], linker_target[0]
    )
    remaining = _dummy_indices(target_linker, atom_map=2)
    if len(remaining) != 1:
        raise ValueError('Expected one remaining linker [*:2] attachment point')
    product = _connect_dummies(target_linker, e3, remaining[0], e3_dummies[0])
    if _dummy_indices(product):
        raise ValueError('Product contains an unconsumed attachment point')
    return Chem.MolToSmiles(product, canonical=True)


def enumerate_linkers(target_fragment_smi, e3_fragment_smi, linkers=None):
    '''Build one connected molecule per valid linker; raise on bad fragments.'''
    linkers = LINKERS if linkers is None else linkers
    return {
        name: build_protac(target_fragment_smi, linker, e3_fragment_smi)
        for name, linker in linkers.items()
    }


def compute_protac_size(protac_smi):
    '''Compute descriptive properties without imposing universal cutoffs.'''
    from rdkit.Chem import Descriptors

    mol = Chem.MolFromSmiles(protac_smi)
    if mol is None:
        raise ValueError('PROTAC SMILES is invalid')
    return {
        'MolWt': Descriptors.MolWt(mol),
        'TPSA': Descriptors.TPSA(mol),
        'LogP': Descriptors.MolLogP(mol),
        'RotBonds': Descriptors.NumRotatableBonds(mol),
    }


if __name__ == '__main__':
    # Demonstration fragments only; each has one explicit exit vector.
    target_fragment = '[*]c1ccc(CN2CCOCC2)cc1'
    e3_fragment = '[*]CC(=O)N1CCC(=O)NC1=O'

    for name, smiles in list(enumerate_linkers(target_fragment, e3_fragment).items())[:5]:
        print(f'{name}: {smiles}')
        print(f'  properties: {compute_protac_size(smiles)}')
