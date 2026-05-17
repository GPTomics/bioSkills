# Reference: RDKit 2024.09+, pandas 2.2+, scikit-learn 1.4+ | Verify API if version differs

from collections import defaultdict
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
import pandas as pd
import random


def get_bemis_murcko_scaffold(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    return Chem.MolToSmiles(scaffold)


def get_generic_framework(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    framework = MurckoScaffold.MakeScaffoldGeneric(scaffold)
    return Chem.MolToSmiles(framework)


def scaffold_clusters(df, smiles_col='smiles'):
    '''Group compound indices by Bemis-Murcko scaffold.'''
    clusters = defaultdict(list)
    for i, row in df.iterrows():
        scaff = get_bemis_murcko_scaffold(row[smiles_col])
        if scaff is not None:
            clusters[scaff].append(i)
    return dict(clusters)


def scaffold_split(df, smiles_col='smiles', train_frac=0.8, seed=42):
    '''Assign whole scaffolds to train or test to prevent leakage in QSAR.'''
    clusters = scaffold_clusters(df, smiles_col)
    # Sort by cluster size desc so large clusters go to train first
    scaffold_sets = sorted(clusters.values(), key=lambda x: len(x), reverse=True)
    n_total = sum(len(s) for s in scaffold_sets)
    n_train = int(n_total * train_frac)

    rng = random.Random(seed)
    rng.shuffle(scaffold_sets[1:])  # keep biggest in train; shuffle rest

    train_idx = []
    test_idx = []
    for scaff_set in scaffold_sets:
        # Prefer adding to train if it fits; else test
        if len(train_idx) + len(scaff_set) <= n_train:
            train_idx.extend(scaff_set)
        else:
            test_idx.extend(scaff_set)
    return df.iloc[train_idx].reset_index(drop=True), df.iloc[test_idx].reset_index(drop=True)


def detect_analog_series(df, smiles_col='smiles', min_size=3):
    '''Identify analog series (scaffold clusters of >= min_size).'''
    clusters = scaffold_clusters(df, smiles_col)
    series = {scaff: cmpds for scaff, cmpds in clusters.items()
              if len(cmpds) >= min_size}
    summary = pd.DataFrame([
        {'scaffold': scaff, 'n_members': len(cmpds)}
        for scaff, cmpds in series.items()
    ]).sort_values('n_members', ascending=False)
    return series, summary


if __name__ == '__main__':
    sample_smiles = [
        'c1ccc(C(=O)NCC)cc1F',
        'c1ccc(C(=O)NCCC)cc1Cl',
        'c1ccc(C(=O)NCC)cc1Br',
        'CCC(=O)N1CCN(c2ccc(F)cc2)CC1',
        'CCC(=O)N1CCN(c2ccc(Cl)cc2)CC1',
    ]
    df = pd.DataFrame({'smiles': sample_smiles, 'pIC50': [6.5, 6.2, 6.0, 7.1, 7.0]})
    train, test = scaffold_split(df, train_frac=0.6, seed=42)
    print(f'Train: {len(train)}, Test: {len(test)}')
    series, summary = detect_analog_series(df, min_size=2)
    print(summary.to_string())
