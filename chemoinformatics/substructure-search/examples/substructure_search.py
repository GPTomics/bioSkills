#!/usr/bin/env python3
'''
Substructure searching with SMARTS patterns.
'''

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D


COMMON_SMARTS = {
    'hydroxyl': '[OH]',
    'primary_amine': '[NH2]',
    'secondary_amine': '[NH1]',
    'carboxylic_acid': '[CX3](=O)[OX2H1]',
    'amide': '[CX3](=O)[NX3]',
    'ester': '[CX3](=O)[OX2][C]',
    'benzene': 'c1ccccc1',
    'aromatic': '[a]',
    'halogen': '[F,Cl,Br,I]',
    'nitro': '[N+]([O-])=O',
    'sulfonamide': '[S](=O)(=O)[NX3]',
    'ketone': '[CX3](=O)[C]',
    'aldehyde': '[CX3H1](=O)',
}


def has_substructure(mol, smarts):
    '''Check if molecule contains substructure.'''
    pattern = Chem.MolFromSmarts(smarts)
    if pattern is None:
        raise ValueError(f'Invalid SMARTS: {smarts}')
    return mol.HasSubstructMatch(pattern)


def get_matches(mol, smarts):
    '''Get all substructure matches as atom indices.'''
    pattern = Chem.MolFromSmarts(smarts)
    if pattern is None:
        raise ValueError(f'Invalid SMARTS: {smarts}')
    return mol.GetSubstructMatches(pattern)


def filter_by_substructure(molecules, smarts, exclude=False):
    '''Filter molecules by substructure presence/absence.'''
    pattern = Chem.MolFromSmarts(smarts)
    if pattern is None:
        raise ValueError(f'Invalid SMARTS: {smarts}')

    filtered = []
    for mol in molecules:
        if mol is None:
            continue
        has_match = mol.HasSubstructMatch(pattern)
        if exclude:
            if not has_match:
                filtered.append(mol)
        else:
            if has_match:
                filtered.append(mol)
    return filtered


def filter_multiple(molecules, include=None, exclude=None):
    '''Filter by multiple inclusion and exclusion patterns.'''
    result = list(molecules)

    if include:
        for smarts in include:
            pattern = Chem.MolFromSmarts(smarts)
            result = [m for m in result if m and m.HasSubstructMatch(pattern)]

    if exclude:
        for smarts in exclude:
            pattern = Chem.MolFromSmarts(smarts)
            result = [m for m in result if m and not m.HasSubstructMatch(pattern)]

    return result


def identify_functional_groups(mol, patterns=None):
    '''Identify functional groups in molecule.'''
    if patterns is None:
        patterns = COMMON_SMARTS

    found = {}
    for name, smarts in patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            found[name] = len(matches)
    return found


def draw_with_highlight(mol, smarts, filename, size=(400, 300)):
    '''Draw molecule with substructure highlighted.'''
    pattern = Chem.MolFromSmarts(smarts)
    match = mol.GetSubstructMatch(pattern)

    drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
    if match:
        drawer.DrawMolecule(mol, highlightAtoms=match)
    else:
        drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    with open(filename, 'wb') as f:
        f.write(drawer.GetDrawingText())


def find_scaffold_matches(molecules, scaffold_smarts):
    '''Find molecules matching a scaffold pattern.'''
    pattern = Chem.MolFromSmarts(scaffold_smarts)
    matches = []
    for i, mol in enumerate(molecules):
        if mol and mol.HasSubstructMatch(pattern):
            matches.append(i)
    return matches


if __name__ == '__main__':
    mol = Chem.MolFromSmiles('c1ccc(O)cc1CC(=O)O')
    print(f'Molecule: {Chem.MolToSmiles(mol)}')

    print('\nFunctional groups found:')
    groups = identify_functional_groups(mol)
    for name, count in groups.items():
        print(f'  {name}: {count}')

    print(f'\nHas hydroxyl: {has_substructure(mol, "[OH]")}')
    print(f'Has primary amine: {has_substructure(mol, "[NH2]")}')

    matches = get_matches(mol, '[OH]')
    print(f'\nHydroxyl positions: {matches}')

    library = [Chem.MolFromSmiles(s) for s in ['CCO', 'CCN', 'CCC', 'c1ccccc1O', 'c1ccccc1N']]
    alcohols = filter_by_substructure(library, '[OH]')
    print(f'\nAlcohols in library: {len(alcohols)}')
