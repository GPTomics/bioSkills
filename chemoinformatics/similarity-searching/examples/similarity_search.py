#!/usr/bin/env python3
'''
Molecular similarity searching and clustering with RDKit.
'''
# Reference: rdkit 2024.09+ | Verify API if version differs

from rdkit import Chem, DataStructs
from rdkit.Chem import MACCSkeys, rdFingerprintGenerator
from rdkit.ML.Cluster import Butina
from rdkit.Chem import rdFMCS


SUPPORTED_FP_TYPES = {'ecfp4', 'maccs'}
ECFP4_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def validate_fp_type(fp_type):
    '''Reject fingerprint names that this example does not implement.'''
    if fp_type not in SUPPORTED_FP_TYPES:
        supported = ', '.join(sorted(SUPPORTED_FP_TYPES))
        raise ValueError(f'Unsupported fingerprint type: {fp_type}. Choose from: {supported}')


def tanimoto_similarity(mol1, mol2, fp_type='ecfp4'):
    '''Calculate Tanimoto similarity between two molecules.'''
    validate_fp_type(fp_type)
    if fp_type == 'ecfp4':
        fp1 = ECFP4_GENERATOR.GetFingerprint(mol1)
        fp2 = ECFP4_GENERATOR.GetFingerprint(mol2)
    elif fp_type == 'maccs':
        fp1 = MACCSkeys.GenMACCSKeys(mol1)
        fp2 = MACCSkeys.GenMACCSKeys(mol2)
    return DataStructs.TanimotoSimilarity(fp1, fp2)


def find_similar(query_smiles, library_smiles, threshold=0.7, fp_type='ecfp4'):
    '''
    Find molecules similar to query in library.

    Returns list of (smiles, similarity) sorted by similarity.
    '''
    validate_fp_type(fp_type)
    query = Chem.MolFromSmiles(query_smiles)
    if query is None:
        raise ValueError('Invalid query SMILES')

    if fp_type == 'ecfp4':
        query_fp = ECFP4_GENERATOR.GetFingerprint(query)
    elif fp_type == 'maccs':
        query_fp = MACCSkeys.GenMACCSKeys(query)

    hits = []
    for smiles in library_smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        if fp_type == 'ecfp4':
            lib_fp = ECFP4_GENERATOR.GetFingerprint(mol)
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
    cutoff = 1 - centroid-neighbor similarity threshold.
    Returned indices refer to the original input sequence.
    '''
    fps = []
    original_indices = []
    for original_idx, mol in enumerate(molecules):
        if mol is not None:
            fps.append(ECFP4_GENERATOR.GetFingerprint(mol))
            original_indices.append(original_idx)

    n = len(fps)
    dists = []
    for i in range(1, n):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1 - s for s in sims])

    clusters = Butina.ClusterData(dists, n, cutoff, isDistData=True)
    return tuple(
        tuple(original_indices[filtered_idx] for filtered_idx in cluster)
        for cluster in clusters
    )


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

    validate_fp_type(fp_type)
    fps = []
    for mol in molecules:
        if fp_type == 'ecfp4':
            fps.append(ECFP4_GENERATOR.GetFingerprint(mol))
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
