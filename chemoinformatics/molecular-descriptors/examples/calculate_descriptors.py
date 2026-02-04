#!/usr/bin/env python3
'''
Molecular descriptor and fingerprint calculation with RDKit.
'''

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem.QED import qed
from rdkit.Chem import MACCSkeys
from rdkit.ML.Descriptors import MoleculeDescriptors
import numpy as np
import pandas as pd


def get_morgan_fingerprint(mol, radius=2, n_bits=2048, use_chirality=False):
    '''
    Generate Morgan fingerprint (ECFP).
    ECFP4 = radius 2, ECFP6 = radius 3
    '''
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits, useChirality=use_chirality)
    return np.array(fp)


def get_maccs_keys(mol):
    '''Generate MACCS keys fingerprint (167 bits).'''
    fp = MACCSkeys.GenMACCSKeys(mol)
    return np.array(fp)


def calculate_lipinski(mol):
    '''Calculate Lipinski Rule of 5 properties.'''
    return {
        'MW': Descriptors.MolWt(mol),
        'LogP': Descriptors.MolLogP(mol),
        'HBD': Lipinski.NumHDonors(mol),
        'HBA': Lipinski.NumHAcceptors(mol)
    }


def passes_lipinski(mol):
    '''Check Lipinski Rule of 5 compliance.'''
    props = calculate_lipinski(mol)
    return (props['MW'] <= 500 and props['LogP'] <= 5 and props['HBD'] <= 5 and props['HBA'] <= 10)


def calculate_druglikeness(mol):
    '''Calculate comprehensive drug-likeness properties.'''
    if mol is None:
        return None

    props = calculate_lipinski(mol)
    props.update({
        'TPSA': Descriptors.TPSA(mol),
        'RotatableBonds': Lipinski.NumRotatableBonds(mol),
        'AromaticRings': Lipinski.NumAromaticRings(mol),
        'QED': qed(mol),
        'LipinskiViolations': sum([
            props['MW'] > 500, props['LogP'] > 5, props['HBD'] > 5, props['HBA'] > 10
        ])
    })

    props['VeberCompliant'] = props['RotatableBonds'] <= 10 and props['TPSA'] <= 140

    return props


def calculate_3d_descriptors(mol):
    '''Calculate 3D descriptors (requires conformer).'''
    from rdkit.Chem import Descriptors3D

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol)

    return {
        'Asphericity': Descriptors3D.Asphericity(mol),
        'Eccentricity': Descriptors3D.Eccentricity(mol),
        'InertialShapeFactor': Descriptors3D.InertialShapeFactor(mol),
        'RadiusOfGyration': Descriptors3D.RadiusOfGyration(mol)
    }


def calculate_all_descriptors(mol):
    '''Calculate all available RDKit descriptors.'''
    descriptor_names = [d[0] for d in Descriptors.descList]
    calculator = MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_names)
    descriptors = calculator.CalcDescriptors(mol)
    return dict(zip(descriptor_names, descriptors))


def batch_calculate_descriptors(molecules, descriptor_list=None):
    '''Calculate descriptors for multiple molecules.'''
    if descriptor_list is None:
        descriptor_list = ['MolWt', 'MolLogP', 'TPSA', 'NumHDonors', 'NumHAcceptors', 'NumRotatableBonds']

    results = []
    for mol in molecules:
        if mol is None:
            results.append({d: None for d in descriptor_list})
            continue
        row = {}
        for name in descriptor_list:
            if name == 'QED':
                row[name] = qed(mol)
            else:
                row[name] = getattr(Descriptors, name)(mol)
        results.append(row)

    return pd.DataFrame(results)


def batch_fingerprints(molecules, fp_type='ecfp4', n_bits=2048):
    '''Generate fingerprints for multiple molecules.'''
    fps = []
    for mol in molecules:
        if mol is None:
            fps.append(np.zeros(n_bits))
            continue
        if fp_type == 'ecfp4':
            fps.append(get_morgan_fingerprint(mol, radius=2, n_bits=n_bits))
        elif fp_type == 'ecfp6':
            fps.append(get_morgan_fingerprint(mol, radius=3, n_bits=n_bits))
        elif fp_type == 'maccs':
            fps.append(get_maccs_keys(mol))
    return np.array(fps)


if __name__ == '__main__':
    mol = Chem.MolFromSmiles('CCO')

    print('Morgan ECFP4:')
    fp = get_morgan_fingerprint(mol)
    print(f'  Shape: {fp.shape}, On bits: {fp.sum()}')

    print('\nLipinski properties:')
    for k, v in calculate_lipinski(mol).items():
        print(f'  {k}: {v}')

    print(f'\nPasses Lipinski: {passes_lipinski(mol)}')

    print('\nDrug-likeness:')
    for k, v in calculate_druglikeness(mol).items():
        print(f'  {k}: {v}')
