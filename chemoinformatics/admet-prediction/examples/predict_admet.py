#!/usr/bin/env python3
'''
ADMET prediction and drug-likeness filtering.
'''
# Reference: rdkit 2024.03+, pandas 2.2+ | Verify API if version differs

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from rdkit.Chem.QED import qed
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams


def load_admetlab_results(csv_path):
    '''Load CSV output produced by the current official ADMETlab 3.0 API workflow.'''
    results = pd.read_csv(csv_path)
    if results.empty:
        raise ValueError('ADMETlab result file contains no predictions')
    return results


def calculate_druglikeness(mol):
    '''Calculate drug-likeness properties.'''
    if mol is None:
        return None

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    tpsa = Descriptors.TPSA(mol)
    rotatable = Lipinski.NumRotatableBonds(mol)

    violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])

    return {
        'MW': mw,
        'LogP': logp,
        'HBD': hbd,
        'HBA': hba,
        'TPSA': tpsa,
        'RotatableBonds': rotatable,
        'QED': qed(mol),
        'LipinskiViolations': violations,
        'VeberCompliant': rotatable <= 10 and tpsa <= 140
    }


def flag_pains(molecules):
    '''
    Flag PAINS patterns without treating them as categorical exclusions.
    Returns (unflagged_molecules, flagged_with_descriptions).
    '''
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)

    clean = []
    flagged = []

    for mol in molecules:
        if mol is None:
            continue
        entry = catalog.GetFirstMatch(mol)
        if entry is None:
            clean.append(mol)
        else:
            flagged.append((mol, entry.GetDescription()))

    return clean, flagged


def flag_structural_alerts(molecules, catalogs=None):
    '''Flag structural alerts using multiple catalogs without discarding molecules.'''
    if catalogs is None:
        catalogs = [
            FilterCatalogParams.FilterCatalogs.PAINS,
            FilterCatalogParams.FilterCatalogs.BRENK,
            FilterCatalogParams.FilterCatalogs.NIH
        ]

    params = FilterCatalogParams()
    for cat in catalogs:
        params.AddCatalog(cat)
    catalog = FilterCatalog(params)

    clean = []
    alerts = []

    for mol in molecules:
        if mol is None:
            continue
        matches = catalog.GetMatches(mol)
        if not matches:
            clean.append(mol)
        else:
            alerts.append((mol, [m.GetDescription() for m in matches]))

    return clean, alerts


def annotate_compounds(molecules):
    '''Annotate compounds; project-specific ranking decisions belong downstream.'''
    results = []

    for mol in molecules:
        if mol is None:
            continue

        props = calculate_druglikeness(mol)
        if props is None:
            continue

        results.append((mol, props))

    return results


def batch_druglikeness(molecules):
    '''Calculate drug-likeness for multiple molecules.'''
    results = []
    for mol in molecules:
        props = calculate_druglikeness(mol)
        if props:
            props['SMILES'] = Chem.MolToSmiles(mol)
            results.append(props)
    return pd.DataFrame(results)


if __name__ == '__main__':
    smiles_list = ['CCO', 'c1ccccc1O', 'CC(=O)Oc1ccccc1C(=O)O']
    molecules = [Chem.MolFromSmiles(s) for s in smiles_list]

    print('Drug-likeness analysis:')
    df = batch_druglikeness(molecules)
    print(df[['SMILES', 'MW', 'LogP', 'QED', 'LipinskiViolations']].to_string())

    print('\nPAINS filtering:')
    unflagged, flagged = flag_pains(molecules)
    print(f'Unflagged: {len(unflagged)}, Flagged: {len(flagged)}')

    print('\nAnnotated compounds:')
    annotated = annotate_compounds(molecules)
    print(f'{len(annotated)} compounds annotated for project-specific review')
