#!/usr/bin/env python3
'''
Molecular similarity searching and clustering with RDKit.
'''

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, MACCSkeys
from rdkit.ML.Cluster import Butina
from rdkit.Chem import rdFMCS


def tanimoto_similarity(mol1, mol2, fp_type='ecfp4'):
    '''Calculate Tanimoto similarity between two molecules.'''
    if fp_type == 'ecfp4':
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
    elif fp_type == 'maccs':
        fp1 = MACCSkeys.GenMACCSKeys(mol1)
        fp2 = MACCSkeys.GenMACCSKeys(mol2)
    return DataStructs.TanimotoSimilarity(fp1, fp2)


def find_similar(query_smiles, library_smiles, threshold=0.7, fp_type='ecfp4'):
    '''
    Find molecules similar to query in library.

    Returns list of (smiles, similarity) sorted by similarity.
    '''
    query = Chem.MolFromSmiles(query_smiles)
    if query is None:
        raise ValueError('Invalid query SMILES')

    if fp_type == 'ecfp4':
        query_fp = AllChem.GetMorganFingerprintAsBitVect(query, 2, nBits=2048)
    elif fp_type == 'maccs':
        query_fp = MACCSkeys.GenMACCSKeys(query)

    hits = []
    for smiles in library_smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        if fp_type == 'ecfp4':
            lib_fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        elif fp_type == 'maccs':
            lib_fp = MACCSkeys.GenMACCSKeys(mol)

        sim = DataStructs.TanimotoSimilarity(query_fp, lib_fp)
        if sim >= threshold:
            hits.append((smiles, sim))

    return sorted(hits, key=lambda x: x[1], reverse=True)


def bulk_similarity_search(query_fp, library_fps, threshold=0.7):
    '''Fast similarity search using bulk operations.'''
    similarities = DataStructs.BulkTanimotoSimilarity(query_fp, library_fps)
    hits = [(i, sim) for i, sim in enumerate(similarities) if sim >= threshold]
    return sorted(hits, key=lambda x: x[1], reverse=True)


def cluster_molecules(molecules, cutoff=0.4):
    '''
    Cluster molecules by Tanimoto similarity using Butina algorithm.
    cutoff = 1 - similarity_threshold (0.4 = 60% similarity threshold)
    '''
    fps = []
    for mol in molecules:
        if mol is not None:
            fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048))

    n = len(fps)
    dists = []
    for i in range(1, n):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1 - s for s in sims])

    clusters = Butina.ClusterData(dists, n, cutoff, isDistData=True)
    return clusters


def find_mcs(molecules, timeout=60):
    '''Find maximum common substructure among molecules.'''
    mcs = rdFMCS.FindMCS(molecules, timeout=timeout, matchValences=False, ringMatchesRingOnly=True)
    return {
        'smarts': mcs.smartsString,
        'num_atoms': mcs.numAtoms,
        'num_bonds': mcs.numBonds
    }


def similarity_matrix(molecules, fp_type='ecfp4'):
    '''Calculate pairwise similarity matrix.'''
    import numpy as np

    fps = []
    for mol in molecules:
        if fp_type == 'ecfp4':
            fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048))
        elif fp_type == 'maccs':
            fps.append(MACCSkeys.GenMACCSKeys(mol))

    n = len(fps)
    sim_matrix = np.ones((n, n))
    for i in range(n):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps)
        sim_matrix[i, :] = sims

    return sim_matrix


if __name__ == '__main__':
    mol1 = Chem.MolFromSmiles('CCO')
    mol2 = Chem.MolFromSmiles('CCCO')
    mol3 = Chem.MolFromSmiles('c1ccccc1')

    print(f'Ethanol vs Propanol: {tanimoto_similarity(mol1, mol2):.3f}')
    print(f'Ethanol vs Benzene: {tanimoto_similarity(mol1, mol3):.3f}')

    library = ['CCO', 'CCCO', 'CCCCO', 'CC(C)O', 'c1ccccc1O']
    hits = find_similar('CCO', library, threshold=0.5)
    print('\nSimilar to ethanol:')
    for smiles, sim in hits:
        print(f'  {smiles}: {sim:.3f}')

    molecules = [Chem.MolFromSmiles(s) for s in library]
    clusters = cluster_molecules(molecules, cutoff=0.3)
    print(f'\nFound {len(clusters)} clusters at 70% similarity')
