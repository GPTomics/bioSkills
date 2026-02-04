#!/usr/bin/env python3
'''
Molecular I/O operations with RDKit.
Reading, writing, and standardizing molecular structures.
'''

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Draw import rdMolDraw2D
from pathlib import Path


def read_sdf(filepath):
    '''Read molecules from SDF file.'''
    supplier = Chem.SDMolSupplier(str(filepath))
    molecules = []
    for mol in supplier:
        if mol is not None:
            molecules.append(mol)
    print(f'Loaded {len(molecules)} molecules from {filepath}')
    return molecules


def read_smiles_file(filepath, delimiter='\t', smiles_col=0, name_col=1):
    '''Read molecules from SMILES file.'''
    molecules = []
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split(delimiter)
            if not parts:
                continue
            smiles = parts[smiles_col]
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                if len(parts) > name_col:
                    mol.SetProp('_Name', parts[name_col])
                molecules.append(mol)
    print(f'Loaded {len(molecules)} molecules from {filepath}')
    return molecules


def write_sdf(molecules, filepath, properties=None):
    '''Write molecules to SDF file.'''
    writer = Chem.SDWriter(str(filepath))
    for mol in molecules:
        if mol is not None:
            writer.write(mol)
    writer.close()
    print(f'Wrote {len(molecules)} molecules to {filepath}')


def write_smiles(molecules, filepath, include_name=True):
    '''Write molecules to SMILES file.'''
    with open(filepath, 'w') as f:
        for mol in molecules:
            if mol is None:
                continue
            smiles = Chem.MolToSmiles(mol)
            if include_name and mol.HasProp('_Name'):
                name = mol.GetProp('_Name')
                f.write(f'{smiles}\t{name}\n')
            else:
                f.write(f'{smiles}\n')


def standardize_molecule(mol):
    '''
    Full standardization pipeline.
    Order: Sanitize -> Normalize -> Neutralize -> Canonicalize tautomer -> Strip salts
    '''
    if mol is None:
        return None

    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None

    normalizer = rdMolStandardize.Normalizer()
    mol = normalizer.normalize(mol)

    uncharger = rdMolStandardize.Uncharger()
    mol = uncharger.uncharge(mol)

    enumerator = rdMolStandardize.TautomerEnumerator()
    mol = enumerator.Canonicalize(mol)

    remover = rdMolStandardize.FragmentRemover()
    mol = remover.remove(mol)

    return mol


def standardize_library(molecules):
    '''Standardize a list of molecules.'''
    standardized = []
    failed = 0
    for mol in molecules:
        std = standardize_molecule(mol)
        if std is not None:
            standardized.append(std)
        else:
            failed += 1
    print(f'Standardized {len(standardized)} molecules, {failed} failed')
    return standardized


def draw_molecule(mol, filename, size=(400, 300)):
    '''Draw molecule to PNG using rdMolDraw2D.'''
    drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    with open(filename, 'wb') as f:
        f.write(drawer.GetDrawingText())


def draw_with_substructure(mol, smarts, filename, size=(400, 300)):
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


def convert_with_openbabel(input_file, output_file, in_format=None, out_format=None):
    '''
    Convert between formats using Open Babel.
    Note: Open Babel 3.x uses 'from openbabel import pybel'
    '''
    from openbabel import pybel

    if in_format is None:
        in_format = Path(input_file).suffix[1:]
    if out_format is None:
        out_format = Path(output_file).suffix[1:]

    output = pybel.Outputfile(out_format, str(output_file), overwrite=True)
    for mol in pybel.readfile(in_format, str(input_file)):
        output.write(mol)
    output.close()


if __name__ == '__main__':
    mol = Chem.MolFromSmiles('CCO')
    print(f'Canonical SMILES: {Chem.MolToSmiles(mol)}')

    mol_with_salt = Chem.MolFromSmiles('CC(=O)O.[Na]')
    standardized = standardize_molecule(mol_with_salt)
    print(f'Standardized: {Chem.MolToSmiles(standardized)}')
