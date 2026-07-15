---
name: bio-scaffold-analysis
description: Analyzes chemical libraries by scaffold using Bemis-Murcko scaffolds, generic frameworks, cyclic skeletons, matched molecular pair (MMP) analysis via mmpdb, R-group decomposition, Free-Wilson analysis, scaffold hopping, and chemotype-aware ML train/test splits. Use when identifying chemotype clusters in a library, deriving SAR transformation rules, decomposing series into R-groups, performing scaffold-balanced QSAR splits, or planning analog campaigns.
tool_type: python
primary_tool: RDKit
---

## Version Compatibility

Reference examples tested with: RDKit 2024.09+, mmpdb 3.1+, scikit-learn 1.4+, datamol 0.12+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Scaffold Analysis

Analyze chemical libraries by their underlying scaffolds. Bemis-Murcko (1996) is the canonical scaffold decomposition: ring systems + linkers, with all R-groups stripped. Generic framework + cyclic skeleton are progressively-more-abstract views. Scaffold analysis underpins QSAR train/test splits (preventing data leakage), library diversity assessment, chemotype clustering, R-group decomposition for SAR modeling, and matched molecular pair analysis (MMPA). The choice of scaffold representation determines whether two compounds are "the same series" -- a critical decision for medicinal chemistry workflows.

For reaction-based enumeration and Free-Wilson, see `chemoinformatics/reaction-enumeration`. For scaffold-hopping via fingerprints, see `chemoinformatics/similarity-searching`. For 3D shape-based scaffold hopping, see `chemoinformatics/shape-similarity`.

## Scaffold Representation Taxonomy

| Representation | Origin | Definition | Use case | Fails when |
|----------------|--------|------------|----------|------------|
| Bemis-Murcko scaffold | Bemis & Murcko 1996 | Ring systems + linkers, R-groups stripped | Default chemotype identifier | Linear molecules (no rings) -> empty scaffold |
| Generic framework | Bemis & Murcko 1996 | Bemis-Murcko with all atoms set to C, all bonds single | Topology comparison | Loses heteroatom info |
| Cyclic skeleton (CSK) | Custom RDKit transformation | Ring atoms only, all C, all single | Pure ring-topology view | Loses linker info; not a built-in Murcko option |
| Murcko atom indices | Derived by matching the scaffold to the parent | Parent-molecule atom indices | Programmatic operations | Symmetry can yield multiple equivalent matches |

```python
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def all_scaffold_views(smi):
    mol = Chem.MolFromSmiles(smi)
    bm = MurckoScaffold.GetScaffoldForMol(mol)
    bm_smi = Chem.MolToSmiles(bm)

    generic = MurckoScaffold.MakeScaffoldGeneric(bm)
    generic_smi = Chem.MolToSmiles(generic)

    return {
        'bemis_murcko': bm_smi,
        'generic_framework': generic_smi,
    }
```

Example: `Cc1ccc(C(=O)NCC2CCCC2)cc1` -> Bemis-Murcko `c1ccc(C(=O)NCC2CCCC2)cc1`; generic `C1CCC(C(C)CCC2CCCC2)CC1` in current RDKit.

## Library Chemotype Clustering

**Goal:** Group compounds by shared Bemis-Murcko scaffold.

**Approach:** Compute scaffold for each compound; group by scaffold SMILES.

```python
from collections import defaultdict

def scaffold_clusters(smiles_list):
    clusters = defaultdict(list)
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        scaffold_smi = Chem.MolToSmiles(scaffold)
        clusters[scaffold_smi].append(smi)
    return clusters
```

Output: dict {scaffold_smiles: [compound_smiles, ...]}. Cluster sizes inform library diversity.

## Bemis-Murcko Scaffold Split (ML)

For QSAR / ML, random train/test split causes data leakage: compounds from the same chemotype (analogs in same series) end up in both. Bemis-Murcko split puts entire scaffolds in train or test, never both.

```python
from rdkit.Chem.Scaffolds import MurckoScaffold

def scaffold_split(df, smiles_col='smiles', train_frac=0.8, seed=42):
    import random
    rng = random.Random(seed)

    scaffolds = defaultdict(list)
    invalid_positions = []
    for pos, smi in enumerate(df[smiles_col].tolist()):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            invalid_positions.append(pos)
            continue
        scaff = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
        scaffolds[scaff].append(pos)

    if invalid_positions:
        raise ValueError(f'Invalid SMILES at row positions: {invalid_positions}')

    scaffold_sets = list(scaffolds.values())
    rng.shuffle(scaffold_sets)
    scaffold_sets.sort(key=lambda x: len(x), reverse=True)

    n_total = sum(len(s) for s in scaffold_sets)
    n_train = int(n_total * train_frac)
    if len(scaffold_sets) < 2:
        raise ValueError('A scaffold split requires at least two scaffolds')
    train_idx = list(scaffold_sets[0])
    test_idx = []
    for i, scaff_set in enumerate(scaffold_sets[1:], start=1):
        if not test_idx and i == len(scaffold_sets) - 1:
            test_idx.extend(scaff_set)
        elif abs(len(train_idx) + len(scaff_set) - n_train) < abs(len(train_idx) - n_train):
            train_idx.extend(scaff_set)
        else:
            test_idx.extend(scaff_set)

    return df.iloc[train_idx], df.iloc[test_idx]
```

**Effect on benchmark metrics:** A scaffold split often produces different performance from a random split because it tests transfer across scaffold groups. The size and meaning of the gap are dataset- and deployment-dependent; it is not a direct universal measure of memorization.

**Caveat:** Bemis-Murcko split is *one* scaffold-split; for production ML, consider time split (newer compounds in test) or activity-cliff-balanced split.

**Class-imbalanced datasets:** Scaffold-only assignment can yield skewed class distributions. Chemprop's `scaffold_balanced` split balances scaffold-group sizes; it is not label-stratified. If both group isolation and label balance are required, use a validated group-aware stratification procedure such as `StratifiedGroupKFold` where its assumptions fit, then audit every fold for scaffold overlap and endpoint balance.

## R-Group Decomposition

**Goal:** Given a defined scaffold and a set of analog compounds, extract the R-group at each numbered attachment point into a tabular SAR matrix.

```python
from rdkit.Chem import rdRGroupDecomposition as rgd

def decompose_series(compounds, scaffold_smiles_with_R):
    scaffold = Chem.MolFromSmiles(scaffold_smiles_with_R)
    if scaffold is None:
        raise ValueError('Invalid scaffold SMARTS/SMILES')
    parsed = [(i, Chem.MolFromSmiles(s)) for i, s in enumerate(compounds)]
    invalid = [i for i, mol in parsed if mol is None]
    if invalid:
        raise ValueError(f'Invalid compound SMILES at positions: {invalid}')
    mols = [mol for _, mol in parsed]
    decomp, unmatched = rgd.RGroupDecompose([scaffold], mols, asSmiles=True)
    unmatched_set = set(unmatched)
    matched_positions = [i for i in range(len(mols)) if i not in unmatched_set]
    return decomp, matched_positions, list(unmatched)

scaffold = 'c1ccc(C(=O)N[*:1])cc1-[*:2]'
compounds = ['c1ccc(C(=O)NCC)cc1F', 'c1ccc(C(=O)NCCC)cc1Cl']
table = decompose_series(compounds, scaffold)
```

Output: list of {'Core': scaffold, 'R1': r1_smiles, 'R2': r2_smiles} dicts. Used for Free-Wilson analysis (see reaction-enumeration skill).

## Matched Molecular Pair Analysis (MMPA) via mmpdb

**Goal:** Mine a SAR dataset for substructure transformations and their associated activity changes.

**Approach:** Fragment all compounds into core + variable side; index pairs differing by one transformation; report delta(activity) per transformation.

```bash
mmpdb fragment data.smi -o data.fragments
mmpdb index data.fragments -o data.mmpdb
mmpdb transform --smiles 'COc1ccccc1' --property pIC50 data.mmpdb
```

Output: ranked transformations with delta(pIC50), N pairs, confidence.

Interpret transformation effects from pair count, chemical-context diversity, dependence among pairs, uncertainty intervals, and prospective validation. Do not convert a universal pair-count/effect-size table into reliability labels.

## Context-Based MMPA

Classical MMPA: "Me -> F always +0.5 log units."
Context-based MMPA: "Me -> F adjacent to amide is +0.5; Me -> F adjacent to ester is -0.1."

Matched-pair effects can depend strongly on the local chemical environment, so report the transformation together with its attachment-point context rather than treating a global mean as universal (Raut & Dixit 2025). Use mmpdb's stored environments or a custom stratified analysis to compare context-specific effects.

## Scaffold Hopping

**Goal:** Find compounds with different scaffold but similar 3D shape / pharmacophore / activity.

| Method | Approach | Tools |
|--------|----------|-------|
| 2D similarity with FCFP4 | Functional-class fingerprint Tanimoto | similarity-searching skill |
| 3D shape (ROCS) | Tanimoto on shape + color volumes | shape-similarity skill |
| Pharmacophore | Common pharmacophore features | pharmacophore-modeling skill |
| Maximum Common Substructure (MCS) | Largest shared substructure | similarity-searching skill (rdFMCS) |
| Deep scaffold hopping | Conditional molecular generation | DeepHop (Zheng et al. 2021) |

For systematic scaffold-hop discovery, combine:
1. Find target's bioactive series
2. Compute 3D pharmacophore from bound conformer
3. ROCS / pharmacophore search against vendor catalogs
4. Filter to compounds with Bemis-Murcko scaffold NOT in training data

## Series Detection

**Goal:** Identify "analog series" within a library -- compounds sharing a scaffold + co-varying R-groups.

```python
def detect_series(smiles_list, min_size=3):
    clusters = scaffold_clusters(smiles_list)
    series = {scaff: cmpds for scaff, cmpds in clusters.items()
              if len(cmpds) >= min_size}
    return series
```

Series counts depend on library provenance, standardization, scaffold definition, and minimum size. Report the observed distribution and use series as one possible unit for SAR analysis.

## Per-Tool Failure Modes

### Bemis-Murcko -- linear molecule yields empty

**Trigger:** Compound has no rings (e.g., fatty acid, simple amine).

**Mechanism:** Bemis-Murcko strips R-groups; no rings = nothing remains.

**Symptom:** Scaffold is empty string; molecules cluster together as "no scaffold".

**Fix:** For linear-rich libraries, augment with linear chain length / functional group features.

### Bemis-Murcko -- spiro / bridged ring confusion

**Trigger:** Compound has spiro or bridged ring system.

**Mechanism:** All ring atoms included; result is the entire ring system without R-groups.

**Symptom:** Apparently different drugs share a "scaffold" because of common spiro center.

**Fix:** Validate visually; use generic framework for topology-only comparison.

### Generic framework -- loses heteroatom info

**Trigger:** Distinguishing pyridine vs benzene scaffolds.

**Mechanism:** `MakeScaffoldGeneric` sets all atoms to C.

**Symptom:** Pyridine and benzene scaffolds reported as identical.

**Fix:** Use Bemis-Murcko (heteroatoms preserved); generic framework for topology only.

### Scaffold split -- imbalanced classes

**Trigger:** Library has many singletons + few large scaffolds.

**Mechanism:** Large scaffolds dominate; greedy assignment puts them in train.

**Symptom:** Test set is mostly singleton scaffolds; metrics misleading.

**Fix:** Use stratified scaffold split (balance test classes); or scaffold-balanced cross-validation.

### MMPA -- low pair count for novel transformations

**Trigger:** Transformation rare in dataset.

**Mechanism:** Need enough pairs to estimate delta(activity).

**Symptom:** Transformation reports N=2 with very large delta.

**Fix:** Report uncertainty and context diversity, avoid overinterpreting sparse transformations, and seek additional matched evidence where appropriate.

### R-group decomposition -- ambiguous match

**Trigger:** Multiple positions in scaffold could match same R-group.

**Mechanism:** Multiple core embeddings, symmetry, and unlabeled attachment choices can yield assignments that differ from the medicinal-chemistry convention.

**Symptom:** R1/R2 columns mixed up.

**Fix:** Specify labeled attachment points, inspect the returned rows and unmatched indices, and use `RGroupDecompositionParameters` for the intended matching/alignment behavior.

## Reconciliation: Scaffold Definition Disagreements

| Concept | Definition A | Definition B | Pick which |
|---------|--------------|--------------|------------|
| Bemis-Murcko scaffold | Atoms in rings + linkers | Same | RDKit default |
| Generic framework | All C, all single bonds | All C, original bonds | `MakeScaffoldGeneric` implements the first; preserve bond orders with an explicit custom transformation |
| Cyclic skeleton | Only ring atoms | Only ring atoms, generic | Implement explicitly; it is not an RDKit Murcko flag |
| "Series" | Same Bemis-Murcko | Tanimoto > 0.8 + same MW | Bemis-Murcko for SAR; Tanimoto for screening |

For ML splits: Bemis-Murcko. For library diversity: Bemis-Murcko + cluster size. For series detection: Bemis-Murcko + R-group decomposition.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Murcko scaffold includes unexpected linker atoms | Bemis-Murcko linkers connect ring systems by definition | Inspect the definition; for hierarchical networks use `rdScaffoldNetwork.ScaffoldNetworkParams` with `CreateScaffoldNetwork` |
| Singleton scaffolds dominate library | Aggressive standardization | Check for tautomer-induced scaffold variation; canonicalize first |
| R-group decomposition empty | Mol doesn't match scaffold | Use FMCS to find actual shared core |
| mmpdb missing transformations | Cores too restrictive | Try smaller core requirement |
| Scaffold split gives all to train | Few scaffolds; large clusters | Add singleton-spread strategy; use Murcko-and-Linker variant |
| Generic framework same for different drugs | Stripped heteroatom info | Use Bemis-Murcko (preserves heteroatoms) |
| MakeScaffoldGeneric error | RDKit version issue | RDKit 2024.09+ uses `Chem.Scaffolds.MurckoScaffold` |

## References

- Bemis GW, Murcko MA. *J. Med. Chem.* 39:2887-2893 (1996) -- original scaffold framework (DOI 10.1021/jm9602928).
- Hu, Stumpfe & Bajorath, *J. Med. Chem.* 60:1238-1246 (2017), DOI 10.1021/acs.jmedchem.6b01437 -- modern scaffold hopping review.
- Hussain J, Rea C. *J. Chem. Inf. Model.* 50:339-348 (2010) -- MMPA core method (DOI 10.1021/ci900450m).
- Raut & Dixit, *RSC Med. Chem.* 16:3281-3290 (2025), DOI 10.1039/D4MD01012D -- local-environment effects in matched molecular pairs.
- Zheng et al., *J. Cheminformatics* 13:87 (2021), DOI 10.1186/s13321-021-00565-5 -- DeepHop conditional scaffold hopping.
- Yang K et al., *J. Chem. Inf. Model.* 59:3370-3388 (2019) -- Chemprop molecular-property prediction (DOI 10.1021/acs.jcim.9b00237).
- RDKit R-group decomposition API: https://www.rdkit.org/docs/source/rdkit.Chem.rdRGroupDecomposition.html
- Chemprop splitting documentation: https://chemprop.readthedocs.io/en/main/tutorial/python/data/splitting.html

## Related Skills

- chemoinformatics/molecular-io - Parse compounds
- chemoinformatics/molecular-standardization - Standardize before scaffold extraction
- chemoinformatics/reaction-enumeration - Free-Wilson analysis on R-decomposition
- chemoinformatics/similarity-searching - 2D scaffold-hopping (FCFP4, AtomPair)
- chemoinformatics/shape-similarity - 3D scaffold-hopping
- chemoinformatics/qsar-modeling - Scaffold-aware splitting for QSAR
- chemoinformatics/generative-design - Scaffold-decoration generative tasks
