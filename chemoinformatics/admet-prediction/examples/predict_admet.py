#!/usr/bin/env python3
'''
ADMET prediction and drug-likeness filtering.
'''

import requests
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from rdkit.Chem.QED import qed
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams


def predict_admetlab(smiles_list, api_url='https://admetlab3.scbdd.com/api/predict'):
    '''
    Predict ADMET using ADMETlab 3.0 API.
    ADMETlab 3.0 provides 119 endpoints with uncertainty.

    Note: SwissADME has NO API (web-only).
    '''
    payload = {'smiles': smiles_list}
    response = requests.post(api_url, json=payload, timeout=60)
    response.raise_for_status()
    return pd.DataFrame(response.json())


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


def filter_pains(molecules):
    '''
    Filter out PAINS (pan-assay interference compounds).
    Returns (clean_molecules, flagged_with_descriptions)
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


def filter_structural_alerts(molecules, catalogs=None):
    '''Filter for structural alerts using multiple catalogs.'''
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


def prioritize_compounds(molecules, max_lipinski_violations=1, min_qed=0.5):
    '''Multi-stage filtering pipeline.'''
    results = []

    for mol in molecules:
        if mol is None:
            continue

        props = calculate_druglikeness(mol)
        if props is None:
            continue

        if props['LipinskiViolations'] > max_lipinski_violations:
            continue

        if not props['VeberCompliant']:
            continue

        if props['QED'] < min_qed:
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
    clean, flagged = filter_pains(molecules)
    print(f'Clean: {len(clean)}, Flagged: {len(flagged)}')

    print('\nPrioritized compounds:')
    prioritized = prioritize_compounds(molecules)
    print(f'{len(prioritized)} compounds passed all filters')
