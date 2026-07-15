---
name: bio-similarity-searching
description: Performs molecular similarity searching using Tanimoto, Tversky, Dice, and cosine coefficients on bit/count fingerprints with explicit choice rules for symmetric vs asymmetric measures, scaffold-hopping vs lead-optimization regimes, activity-cliff diagnosis, and large-library nearest-neighbor methods (BulkTanimoto, MHFP6 LSH forest, USRCAT). Use when ranking compounds by structural resemblance to a query, clustering libraries, finding analogs, or diagnosing activity cliffs.
tool_type: python
primary_tool: RDKit
---

## Version Compatibility

Reference examples tested with: RDKit 2024.09+, scikit-learn 1.4+, mhfp 1.9+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Similarity Searching

Find structurally similar compounds and cluster libraries by similarity. The choice of similarity coefficient and fingerprint is **task-aware**: Tanimoto for symmetric similarity in lead optimization, Tversky for asymmetric "substructure-like" queries, Dice for higher sensitivity in low-similarity regimes, and MaxCommon Substructure (MCS) for scaffold-hopping. Tanimoto similarity above 0.7 is not a guarantee of activity preservation; activity cliffs (similar molecules with dissimilar activities) are common (Maggiora 2014).

For fingerprint choice, see `chemoinformatics/molecular-descriptors`. For 3D shape similarity, see `chemoinformatics/shape-similarity`.

## Similarity Coefficient Taxonomy

| Coefficient | Formula | Range | Symmetric | Use case | Fails when |
|-------------|---------|-------|-----------|----------|------------|
| Tanimoto | c / (a + b - c) | 0-1 | Yes | Default for ECFP4 similarity, ranking analogs | Saturates at low similarity (drug vs natural product) |
| Dice | 2c / (a + b) | 0-1 | Yes | Bit or nonnegative sparse-count vectors when Dice semantics are intended | Thresholds depend on vector type; analog choice subjective |
| Cosine (Ochiai) | c / sqrt(a*b) | 0-1 | Yes | Count vectors, weighted similarity | Not standard for bit vectors |
| Tversky alpha,beta | c / (alpha*(a-c) + beta*(b-c) + c) | 0-1 | No when alpha != beta | Asymmetric "is A a substructure of B" queries | Parameter choice subjective; alpha=1,beta=0 = substructure-like |
| Hamming | (a + b - 2c) / nBits | 0-1 | Yes | Binary bit vectors when bit-wise disagreement matters | Does not preserve count magnitude |
| Russell-Rao | c / nBits | 0-1 | Yes | Sparse fingerprints | Biased by fingerprint density |
| Kulczynski | (c/a + c/b) / 2 | 0-1 | Yes | When fingerprints have very different bit-density | Less standard |

Where a = set bits in fp1, b = set bits in fp2, c = bits in common.

## When to Use Which Coefficient

| Scenario | Coefficient | Why |
|----------|-------------|-----|
| Standard analog search (drug-like, ECFP4) | Tanimoto, start near 0.7 | Repository starting heuristic; calibrate against project actives and analog judgments |
| Sensitive search at lower similarity | Dice, threshold 0.45 | Dice is roughly 2*Tanimoto/(1+Tanimoto); more sensitive in middle range |
| Substructure-like ranking | Tversky alpha=1, beta=0 | Asymmetric: rewards compounds containing query features |
| Count fingerprints (neural, atom-environment) | Cosine | Bit-vector Tanimoto loses information |
| Activity-cliff diagnosis | Tanimoto + property difference | Detect ECFP4>=0.85 but |delta(activity)|>=2 log units |
| Cross-target / scaffold-hopping | FCFP4 Tanimoto OR AtomPair Tanimoto | Pharmacophore-equivalent matches different scaffolds |
| Metabolomics / natural products | MHFP6 Jaccard | ECFP4 saturates near 0.2 across diverse classes |
| 3D shape | Tanimoto on shape volume overlap | See shape-similarity skill |

## Tanimoto Thresholds (Repository Starting Heuristics)

| Threshold | Interpretation | Caveat |
|-----------|----------------|--------|
| >=0.85 | Likely same scaffold + close analog | Activity cliffs still possible |
| 0.70-0.85 | Same series, R-group variation | Standard "similar" threshold |
| 0.55-0.70 | Related chemotype, different decoration | Useful for series expansion |
| 0.35-0.55 | Distant analog, possible scaffold hop | Many false positives |
| <0.35 | Mostly noise; use 3D shape or pharmacophore instead | ECFP4 not informative |

These bands are working defaults for ECFP4-like fingerprints, not transferable calibration. Inspect the target dataset's similarity distribution and known series before setting a cutoff. Maggiora's similarity principle states "similar molecules tend to have similar activity" -- but **activity cliffs** (Stumpfe & Bajorath 2012) violate this. Treat high ECFP4 similarity as a prioritization signal, not evidence that activity will be preserved.

## Decision Tree by Scenario

| Goal | Workflow | Tools |
|------|----------|-------|
| Find analogs of a hit (lead opt) | ECFP4 Tanimoto >=0.7 search | RDKit `BulkTanimotoSimilarity` |
| Find scaffold hops | FCFP4 OR AtomPair Tanimoto >=0.5 + filter MCS | RDKit + rdFMCS |
| Cluster library by chemotype | Butina clustering at Tanimoto 0.6 cutoff | RDKit `Butina.ClusterData` |
| Diversity sampling | MaxMin selection on Tanimoto | RDKit `rdSimDivPickers.MaxMinPicker` |
| Nearest neighbors in >1M library | LSH forest with MHFP6 | `mhfp.lsh_forest.LSHForestHelper` |
| Activity cliff diagnosis | Tanimoto + pIC50 delta scatter | Custom analysis |
| 3D similarity (shape) | USRCAT / Open3DAlign / ROCS | shape-similarity skill |

## Tanimoto Similarity (single query, large library)

**Goal:** Rank a library by ECFP4 Tanimoto similarity to a query molecule, returning hits above a threshold.

**Approach:** Generate ECFP4 fingerprints for all molecules once, then use `BulkTanimotoSimilarity` for O(N) lookup.

```python
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

def precompute_fps(smiles_list, radius=2, nBits=2048):
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=radius, fpSize=nBits)
    fps = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            fps.append(None)
        else:
            fps.append(generator.GetFingerprint(mol))
    return fps

def search(query_smi, library_fps, threshold=0.7):
    qmol = Chem.MolFromSmiles(query_smi)
    if qmol is None:
        raise ValueError('invalid query SMILES')
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    qfp = generator.GetFingerprint(qmol)
    valid = [(i, fp) for i, fp in enumerate(library_fps) if fp is not None]
    sims = DataStructs.BulkTanimotoSimilarity(qfp, [fp for _, fp in valid])
    return [(source_i, sim) for (source_i, _), sim in zip(valid, sims)
            if sim >= threshold]
```

## Tversky for Asymmetric Substructure-Like Search

**Goal:** Rank a library by how much each compound "contains" the features of the query (asymmetric).

**Approach:** Tversky with alpha=1, beta=0 rewards compounds containing query bits (substructure-like) while ignoring extra bits in the compound.

```python
from rdkit import DataStructs

def tversky_substructure_like(qfp, lib_fps, alpha=1.0, beta=0.0):
    return [DataStructs.TverskySimilarity(qfp, f, alpha, beta) for f in lib_fps if f]
```

Use case: identifying analogs that extend a pharmacophore vs. exact-similarity ranking.

## Butina Clustering

**Goal:** Group a library around Taylor-Butina centroids whose assigned neighbors are within the selected distance cutoff.

**Approach:** Compute upper-triangle distance matrix, apply Taylor-Butina with chosen distance cutoff.

```python
from rdkit.ML.Cluster import Butina

def cluster(mols, cutoff=0.4):
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = [generator.GetFingerprint(m) for m in mols]
    n = len(fps)
    dists = []
    for i in range(1, n):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1 - s for s in sims])
    return Butina.ClusterData(dists, n, cutoff, isDistData=True)
```

`cutoff=0.4` means each assigned member was a neighbor of its selected centroid at Tanimoto >= 0.6. It does **not** guarantee that every pair of non-centroid members has Tanimoto >= 0.6. The first molecule in each returned cluster is the cluster centroid.

**Trade-off:** Butina materializes O(N^2) pairwise distances. Benchmark memory and runtime on the actual library; for much larger collections, use an approximate method such as an MHFP6 LSH forest.

## Diversity Selection (MaxMin)

**Goal:** Select N diverse compounds from a library by maximizing the minimum pairwise distance.

```python
from rdkit.SimDivFilters import rdSimDivPickers

picker = rdSimDivPickers.MaxMinPicker()
n_pick = 100
n_lib = len(fps)
selected = picker.LazyBitVectorPick(fps, n_lib, n_pick, seed=42)
```

`LazyBitVectorPick` is memory-efficient (does not materialize full distance matrix).

## Maximum Common Substructure

**Goal:** Find the largest substructure shared across a set of molecules.

**Approach:** `rdFMCS.FindMCS` with parameters controlling atom/bond equivalence.

```python
from rdkit.Chem import rdFMCS

def mcs_smarts(mols, timeout=60, ring_match='strict', atom_match='elements'):
    params = rdFMCS.MCSParameters()
    params.Timeout = timeout
    if ring_match == 'strict':
        params.BondCompareParameters.MatchFusedRings = True
        params.BondCompareParameters.MatchFusedRingsStrict = True
        params.BondCompareParameters.RingMatchesRingOnly = True
    if atom_match == 'elements':
        params.AtomCompareParameters.MatchValences = False
    result = rdFMCS.FindMCS(mols, params)
    return result.smartsString, result.numAtoms, result.numBonds
```

Use cases: identify scaffold across a series, build scaffold hopping queries, generate consensus pharmacophore.

**Limit:** MCS search can become combinatorial as input count, size, and structural divergence grow. Set a finite timeout, inspect `result.canceled`, and consider pre-clustering or reducing the comparison set; no molecule-count or atom-count boundary guarantees tractability.

## Activity Cliff Diagnosis

**Goal:** Detect pairs of similar molecules with dissimilar activities (cliffs).

**Approach:** Compute pairwise ECFP4 Tanimoto + pIC50 delta. Flag pairs with high similarity and large activity gap.

```python
def activity_cliffs(df, sim_threshold=0.85, activity_gap=2.0, activity_col='pIC50'):
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    mols = [Chem.MolFromSmiles(s) for s in df['smiles']]
    if any(mol is None for mol in mols):
        raise ValueError('activity-cliff input contains invalid SMILES')
    fps = [generator.GetFingerprint(mol) for mol in mols]
    cliffs = []
    for i in range(len(fps)):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[i+1:])
        for j_off, sim in enumerate(sims):
            j = i + 1 + j_off
            if sim >= sim_threshold:
                gap = abs(df[activity_col].iloc[i] - df[activity_col].iloc[j])
                if gap >= activity_gap:
                    cliffs.append((i, j, sim, gap))
    return cliffs
```

Activity cliffs flag (a) measurement noise, (b) cryptic SAR (e.g. ring-flip changing dihedral), (c) protein conformational selection, or (d) actually informative SAR. Cliffs are an opportunity for medchem investigation, not necessarily an error.

## Large-Library Nearest Neighbor (MHFP6 + LSH Forest)

For large libraries, direct all-pairs comparison becomes expensive. The `mhfp` package provides an LSH-forest helper for approximate nearest-neighbor retrieval over MHFP6 fingerprints. Measure recall and latency against an exact subset for the project dataset.

```python
from mhfp.encoder import MHFPEncoder
from mhfp.lsh_forest import LSHForestHelper

encoder = MHFPEncoder(2048)

def build_index(smiles_list):
    forest = LSHForestHelper()
    fingerprints = []
    for i, smiles in enumerate(smiles_list):
        fp = encoder.encode(smiles, radius=3)
        fingerprints.append(fp)
        forest.add(i, fp)
    forest.index()
    return forest, fingerprints

def query_index(forest, qmol, fingerprints, k=10):
    qfp = encoder.encode_mol(qmol, radius=3)
    return forest.query(qfp, k=k, data=fingerprints)
```

The returned neighbors are approximate in MHFP6 space. Benchmark them against an exact MHFP-distance search on a representative subset before choosing LSH parameters.

## Per-Tool Failure Modes

### ECFP4 Tanimoto -- saturation in diverse libraries

**Trigger:** Library spans drug-like + natural products + peptides + metabolites.

**Mechanism:** A local circular fingerprint may not preserve the distinctions needed for a particular mixed-modality retrieval task.

**Symptom:** Known relevant neighbors are not enriched above background, or retrieval metrics and neighborhood stability are poor on target-relevant controls. A low mean pairwise similarity alone is not a universal saturation test.

**Fix:** Benchmark ECFP4 against alternatives such as MHFP6 or MAP4 using held-out analog recovery, scaffold-aware retrieval, or another task-aligned metric. Do not transfer raw-score thresholds between fingerprint families.

### Butina clustering -- O(N^2) memory blowup

**Trigger:** Library >100k molecules.

**Mechanism:** Butina requires upper-triangle distance matrix, ~5e9 floats for 100k compounds.

**Symptom:** OOM error or hours of CPU.

**Fix:** Use approximate clustering (HDBSCAN on UMAP-reduced fingerprints) or LSH-based clustering on MHFP6.

### MCS -- exponential timeout

**Trigger:** Mol set with low overlap, large molecules, or many input mols.

**Mechanism:** MCS search is NP-hard; algorithm tries every atom-mapping permutation within timeout.

**Symptom:** Returns small partial MCS or empty result.

**Fix:** Raise `timeout`; reduce input mol count; pre-cluster by Tanimoto first then MCS within clusters.

### Tanimoto = 1.0 != same molecule

**Trigger:** Comparing fingerprints between two molecules that hash to the same bits.

**Mechanism:** A folded hashed fingerprint can map distinct atom environments to the same bits; collision frequency depends on molecule size, radius, and fingerprint length.

**Symptom:** Two structurally different molecules report Tanimoto 1.0.

**Fix:** For exact identity, compare canonical SMILES or InChIKey, not fingerprint. Use unhashed sparse fingerprint to disambiguate.

### Similarity threshold transfer fails

**Trigger:** Threshold tuned on ECFP4 applied to RDKit FP, AtomPair, or MACCS.

**Mechanism:** Bit-density and fragment-resolution differ; Tanimoto distributions shift.

**Symptom:** "Similar" set is much larger or smaller than expected.

**Fix:** Re-tune the threshold per fingerprint and dataset. AtomPair ~0.55, MACCS ~0.85, ECFP4 ~0.7, and FCFP4 ~0.6 are repository starting heuristics, not universal equivalents.

## Reconciliation: Cliffs Across Methods

If a pair flags as an activity cliff under one representation but not another, treat that as representation sensitivity. Inspect atom mappings, fingerprint environments, assay uncertainty, and the exact structural change; disagreement alone does not establish which substituent caused the activity difference.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `BulkTanimotoSimilarity` output is treated as bit counts | The API returns similarity values, normally floats | Keep the returned values as similarities; inspect input vector types if the output is unexpected |
| Reported similarity > 1 | Custom formula, malformed data, negative features, or an incorrectly normalized external implementation | Verify the coefficient definition and inputs; standard nonnegative RDKit Tanimoto and Tversky similarities are bounded by 1 |
| Cluster centroids change when input order changes | Taylor-Butina tie handling and assignment are order-sensitive | Standardize and sort inputs by a stable identifier before clustering; record `reordering` and the input order |
| MaxMinPicker returns first N inputs | All-zero initial similarity matrix | Seed picker explicitly: `picker.LazyBitVectorPick(fps, n_lib, n_pick, seed=42)` |
| Activity cliff "false positives" | Bit-collisions inflate similarity | Use sparse Morgan or compare canonical SMILES for exact ID |
| Diverse subset has duplicates | Standardization not applied | Canonicalize via `chemoinformatics/molecular-standardization` first |
| Tanimoto incompatible with neural fingerprint | Continuous-valued fingerprint | Use cosine or sklearn `cdist` with `'cosine'` metric |

## References

- Bajorath, *Nat. Rev. Drug Discov.* 1:882-894 (2002) -- integration of virtual and high-throughput screening. https://doi.org/10.1038/nrd941
- Maggiora et al., *J. Med. Chem.* 57:3186-3204 (2014) -- molecular similarity in drug discovery. https://doi.org/10.1021/jm401411z
- Stumpfe & Bajorath, *J. Med. Chem.* 55:2932-2942 (2012) -- activity cliffs. https://doi.org/10.1021/jm201706b
- Tversky, *Psychol. Rev.* 84:327-352 (1977) -- features of similarity (Tversky coefficient). https://doi.org/10.1037/0033-295X.84.4.327
- Probst & Reymond, *J. Cheminformatics* 10:66 (2018) -- MHFP6 MinHash fingerprint. https://doi.org/10.1186/s13321-018-0321-8
- O'Boyle & Sayle, *J. Cheminformatics* 8:36 (2016) -- fingerprint similarity benchmark. https://doi.org/10.1186/s13321-016-0148-0
- Butina, *J. Chem. Inf. Comput. Sci.* 39:747-750 (1999) -- Taylor-Butina clustering. https://doi.org/10.1021/ci9803381
- RDKit, `rdFMCS` API documentation -- MCS ring-comparison parameter names and semantics. https://www.rdkit.org/docs/source/rdkit.Chem.rdFMCS.html

## Related Skills

- chemoinformatics/molecular-descriptors - Generate fingerprints for similarity
- chemoinformatics/molecular-standardization - Canonicalize before comparing
- chemoinformatics/substructure-search - SMARTS pattern-based searching
- chemoinformatics/scaffold-analysis - Scaffold-based similarity (Bemis-Murcko, MMPA)
- chemoinformatics/shape-similarity - 3D shape similarity (USRCAT, ROCS)
- machine-learning/biomarker-discovery - ML on similarity features
