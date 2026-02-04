---
name: bio-similarity-searching
description: Performs molecular similarity searches using Tanimoto coefficient on fingerprints via RDKit. Finds structurally similar compounds using ECFP or MACCS keys and clusters molecules by structural similarity using Butina clustering. Use when finding analogs of a query compound or clustering chemical libraries.
tool_type: python
primary_tool: RDKit
---

# Similarity Searching

Find structurally similar molecules and cluster compound libraries.

## Tanimoto Similarity

```python
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

# Generate fingerprints
mol1 = Chem.MolFromSmiles('CCO')
mol2 = Chem.MolFromSmiles('CCCO')

fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius=2, nBits=2048)
fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, radius=2, nBits=2048)

# Tanimoto similarity (0-1)
similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
print(f'Tanimoto similarity: {similarity:.3f}')
```

## Similarity Thresholds

| Threshold | Interpretation |
|-----------|----------------|
| > 0.85 | Very similar (likely same scaffold) |
| > 0.70 | Similar (likely related series) |
| > 0.50 | Moderate similarity |
| < 0.50 | Dissimilar |

## Search Library Against Query

```python
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

def find_similar_molecules(query_smiles, library, threshold=0.7, fp_type='ecfp4'):
    '''
    Find molecules similar to query in library.

    Args:
        query_smiles: Query molecule SMILES
        library: List of (smiles, name) tuples or SMILES list
        threshold: Minimum Tanimoto similarity
        fp_type: 'ecfp4', 'ecfp6', or 'maccs'
    '''
    query = Chem.MolFromSmiles(query_smiles)
    if query is None:
        raise ValueError('Invalid query SMILES')

    # Generate query fingerprint
    if fp_type == 'ecfp4':
        query_fp = AllChem.GetMorganFingerprintAsBitVect(query, 2, nBits=2048)
    elif fp_type == 'ecfp6':
        query_fp = AllChem.GetMorganFingerprintAsBitVect(query, 3, nBits=2048)
    else:  # maccs
        from rdkit.Chem import MACCSkeys
        query_fp = MACCSkeys.GenMACCSKeys(query)

    # Search library
    hits = []
    for item in library:
        smiles = item[0] if isinstance(item, tuple) else item
        name = item[1] if isinstance(item, tuple) and len(item) > 1 else smiles

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue

        if fp_type == 'ecfp4':
            lib_fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        elif fp_type == 'ecfp6':
            lib_fp = AllChem.GetMorganFingerprintAsBitVect(mol, 3, nBits=2048)
        else:
            lib_fp = MACCSkeys.GenMACCSKeys(mol)

        sim = DataStructs.TanimotoSimilarity(query_fp, lib_fp)
        if sim >= threshold:
            hits.append((smiles, name, sim))

    return sorted(hits, key=lambda x: x[2], reverse=True)
```

## Bulk Similarity Search

```python
from rdkit import DataStructs

def bulk_similarity_search(query_fp, library_fps, threshold=0.7):
    '''
    Fast similarity search using bulk operations.

    Args:
        query_fp: Query fingerprint
        library_fps: List of library fingerprints
        threshold: Minimum similarity
    '''
    # BulkTanimotoSimilarity is faster for large libraries
    similarities = DataStructs.BulkTanimotoSimilarity(query_fp, library_fps)

    hits = [(i, sim) for i, sim in enumerate(similarities) if sim >= threshold]
    return sorted(hits, key=lambda x: x[1], reverse=True)
```

## Butina Clustering

```python
from rdkit import Chem
from rdkit.ML.Cluster import Butina

def cluster_molecules(molecules, cutoff=0.4):
    '''
    Cluster molecules by Tanimoto similarity using Taylor-Butina algorithm.

    Args:
        molecules: List of RDKit mol objects
        cutoff: Distance cutoff (1 - similarity threshold)
               cutoff=0.4 means similarity threshold of 0.6
    '''
    # Generate fingerprints
    fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)
           for m in molecules if m is not None]

    # Calculate distance matrix (upper triangle)
    n = len(fps)
    dists = []
    for i in range(1, n):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1 - s for s in sims])

    # Cluster
    clusters = Butina.ClusterData(dists, n, cutoff, isDistData=True)

    return clusters

# Usage
# clusters = cluster_molecules(molecules, cutoff=0.3)  # 70% similarity
# print(f'Found {len(clusters)} clusters')
# for i, cluster in enumerate(clusters[:5]):
#     print(f'Cluster {i}: {len(cluster)} molecules')
```

## Maximum Common Substructure

```python
from rdkit.Chem import rdFMCS

def find_mcs(molecules, timeout=60):
    '''Find maximum common substructure.'''
    mcs = rdFMCS.FindMCS(
        molecules,
        timeout=timeout,
        matchValences=False,
        ringMatchesRingOnly=True
    )
    return mcs.smartsString, mcs.numAtoms, mcs.numBonds

# Get MCS as molecule for visualization
mcs_smarts, n_atoms, n_bonds = find_mcs(molecules)
mcs_mol = Chem.MolFromSmarts(mcs_smarts)
```

## Related Skills

- molecular-descriptors - Generate fingerprints for similarity
- substructure-search - Pattern-based searching
- molecular-io - Load molecules for searching
