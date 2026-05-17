---
name: bio-similarity-searching
description: Performs molecular similarity searching using Tanimoto, Tversky, Dice, and cosine coefficients on bit/count fingerprints with explicit choice rules for symmetric vs asymmetric measures, scaffold-hopping vs lead-optimization regimes, activity-cliff diagnosis, and large-library nearest-neighbor methods (BulkTanimoto, Annoy MHFP6, USRCAT). Use when ranking compounds by structural resemblance to a query, clustering libraries, finding analogs, or diagnosing activity cliffs.
tool_type: python
primary_tool: RDKit
---

## Version Compatibility

Reference examples tested with: RDKit 2024.09+, scikit-learn 1.4+, annoy 1.17+, mhfp 1.9+.

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
| Dice | 2c / (a + b) | 0-1 | Yes | More sensitive than Tanimoto in 0.3-0.5 range | Bit-vector only; analog choice subjective |
| Cosine (Ochiai) | c / sqrt(a*b) | 0-1 | Yes | Count vectors, weighted similarity | Not standard for bit vectors |
| Tversky alpha,beta | c / (alpha*(a-c) + beta*(b-c) + c) | 0-1 | No when alpha != beta | Asymmetric "is A a substructure of B" queries | Parameter choice subjective; alpha=1,beta=0 = substructure-like |
| Hamming | (a + b - 2c) / nBits | 0-1 | Yes | Count vectors, when bit-wise distance matters | Bit-vector only loses scale |
| Russell-Rao | c / nBits | 0-1 | Yes | Sparse fingerprints | Biased by fingerprint density |
| Kulczynski | (c/a + c/b) / 2 | 0-1 | Yes | When fingerprints have very different bit-density | Less standard |

Where a = set bits in fp1, b = set bits in fp2, c = bits in common.

## When to Use Which Coefficient

| Scenario | Coefficient | Why |
|----------|-------------|-----|
| Standard analog search (drug-like, ECFP4) | Tanimoto, threshold 0.7 | Industry default; calibrated against medchem judgment |
| Sensitive search at lower similarity | Dice, threshold 0.45 | Dice is roughly 2*Tanimoto/(1+Tanimoto); more sensitive in middle range |
| Substructure-like ranking | Tversky alpha=1, beta=0 | Asymmetric: rewards compounds containing query features |
| Count fingerprints (neural, atom-environment) | Cosine | Bit-vector Tanimoto loses information |
| Activity-cliff diagnosis | Tanimoto + property difference | Detect ECFP4>=0.85 but |delta(activity)|>=2 log units |
| Cross-target / scaffold-hopping | FCFP4 Tanimoto OR AtomPair Tanimoto | Pharmacophore-equivalent matches different scaffolds |
| Metabolomics / natural products | MHFP6 Jaccard | ECFP4 saturates near 0.2 across diverse classes |
| 3D shape | Tanimoto on shape volume overlap | See shape-similarity skill |

## Tanimoto Thresholds (calibrated against medchem judgment)

| Threshold | Interpretation | Caveat |
|-----------|----------------|--------|
| >=0.85 | Likely same scaffold + close analog | Activity cliffs still possible |
| 0.70-0.85 | Same series, R-group variation | Standard "similar" threshold |
| 0.55-0.70 | Related chemotype, different decoration | Useful for series expansion |
| 0.35-0.55 | Distant analog, possible scaffold hop | Many false positives |
| <0.35 | Mostly noise; use 3D shape or pharmacophore instead | ECFP4 not informative |

Maggiora's similarity principle states "similar molecules tend to have similar activity" -- but **activity cliffs** (Stumpfe & Bajorath 2012) violate this. Roughly 10-20% of activity-cliff pairs have ECFP4 Tanimoto >=0.7 with delta(pIC50) >=2.

## Decision Tree by Scenario

| Goal | Workflow | Tools |
|------|----------|-------|
| Find analogs of a hit (lead opt) | ECFP4 Tanimoto >=0.7 search | RDKit `BulkTanimotoSimilarity` |
| Find scaffold hops | FCFP4 OR AtomPair Tanimoto >=0.5 + filter MCS | RDKit + rdFMCS |
| Cluster library by chemotype | Butina clustering at Tanimoto 0.6 cutoff | RDKit `Butina.ClusterData` |
| Diversity sampling | MaxMin selection on Tanimoto | RDKit `rdSimDivPickers.MaxMinPicker` |
| Nearest neighbors in >1M library | LSH (MinHash) with MHFP6 | mhfp + Annoy |
| Activity cliff diagnosis | Tanimoto + pIC50 delta scatter | Custom analysis |
| 3D similarity (shape) | USRCAT / Open3DAlign / ROCS | shape-similarity skill |

## Tanimoto Similarity (single query, large library)

**Goal:** Rank a library by ECFP4 Tanimoto similarity to a query molecule, returning hits above a threshold.

**Approach:** Generate ECFP4 fingerprints for all molecules once, then use `BulkTanimotoSimilarity` for O(N) lookup.

```python
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

def precompute_fps(smiles_list, radius=2, nBits=2048):
    fps = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            fps.append(None)
        else:
            fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits))
    return fps

def search(query_smi, library_fps, threshold=0.7):
    qmol = Chem.MolFromSmiles(query_smi)
    qfp = AllChem.GetMorganFingerprintAsBitVect(qmol, 2, nBits=2048)
    sims = DataStructs.BulkTanimotoSimilarity(qfp, [f for f in library_fps if f])
    return [(i, s) for i, s in enumerate(sims) if s >= threshold]
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

**Goal:** Group a library into clusters where intra-cluster Tanimoto >= 1 - cutoff.

**Approach:** Compute upper-triangle distance matrix, apply Taylor-Butina with chosen distance cutoff.

```python
from rdkit.ML.Cluster import Butina

def cluster(mols, cutoff=0.4):
    fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) for m in mols]
    n = len(fps)
    dists = []
    for i in range(1, n):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1 - s for s in sims])
    return Butina.ClusterData(dists, n, cutoff, isDistData=True)
```

`cutoff=0.4` means clusters share Tanimoto >= 0.6. The first molecule in each returned cluster is the cluster centroid.

**Trade-off:** Butina is O(N^2) and scales to ~100k molecules. For 1M+ libraries, use approximate nearest-neighbor (Annoy with MHFP6).

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
        params.BondCompareParameters.MatchRingFusion = True
        params.BondCompareParameters.RingMatchesRingOnly = True
    if atom_match == 'elements':
        params.AtomCompareParameters.MatchValences = False
    result = rdFMCS.FindMCS(mols, params)
    return result.smartsString, result.numAtoms, result.numBonds
```

Use cases: identify scaffold across a series, build scaffold hopping queries, generate consensus pharmacophore.

**Limit:** MCS scales exponentially with molecule overlap. For >50 molecules or molecules >30 atoms, raise `timeout` and accept partial results.

## Activity Cliff Diagnosis

**Goal:** Detect pairs of similar molecules with dissimilar activities (cliffs).

**Approach:** Compute pairwise ECFP4 Tanimoto + pIC50 delta. Flag pairs with high similarity and large activity gap.

```python
def activity_cliffs(df, sim_threshold=0.85, activity_gap=2.0, activity_col='pIC50'):
    fps = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), 2, nBits=2048)
           for s in df['smiles']]
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

## Large-Library Nearest Neighbor (MHFP6 + Annoy)

For libraries >1M compounds, direct pairwise Tanimoto is impractical. MHFP6 (MinHashed fingerprint) + Annoy (LSH) gives approximate top-k in seconds.

```python
from mhfp.encoder import MHFPEncoder
from annoy import AnnoyIndex

encoder = MHFPEncoder(2048)

def build_index(mols, fname):
    idx = AnnoyIndex(2048, 'hamming')
    for i, mol in enumerate(mols):
        fp = encoder.encode(mol, radius=3)
        idx.add_item(i, fp)
    idx.build(50)
    idx.save(fname)
    return idx

def query_index(idx, qmol, k=10):
    qfp = encoder.encode(qmol, radius=3)
    return idx.get_nns_by_vector(qfp, k, include_distances=True)
```

Annoy distance: Hamming on MinHash, related to Jaccard. Trade-off: ~95% recall on top-10 nearest neighbors vs full-Tanimoto search.

## Per-Tool Failure Modes

### ECFP4 Tanimoto -- saturation in diverse libraries

**Trigger:** Library spans drug-like + natural products + peptides + metabolites.

**Mechanism:** ECFP4 bits sparsely cover macrocycles, peptides, glycans; most pairwise Tanimoto reports 0.1-0.3.

**Symptom:** Pairwise Tanimoto distribution narrow; ranking is dominated by noise. Quantitative diagnostic: if mean pairwise Tanimoto on a random sample of N>=1000 from the library is < 0.20, ECFP4 is saturated.

**Fix:** Use MHFP6 (MinHash) or MAP4 (atom-pair MinHash); calibrated for chemical diversity. Re-check mean Tanimoto > 0.30 with new fingerprint.

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

**Mechanism:** Hashed Morgan with nBits=2048 has bit-collision rate ~1-5% for drug-like molecules.

**Symptom:** Two structurally different molecules report Tanimoto 1.0.

**Fix:** For exact identity, compare canonical SMILES or InChIKey, not fingerprint. Use unhashed sparse fingerprint to disambiguate.

### Similarity threshold transfer fails

**Trigger:** Threshold tuned on ECFP4 applied to RDKit FP, AtomPair, or MACCS.

**Mechanism:** Bit-density and fragment-resolution differ; Tanimoto distributions shift.

**Symptom:** "Similar" set is much larger or smaller than expected.

**Fix:** Re-tune threshold per fingerprint. AtomPair similar threshold ~0.55; MACCS ~0.85; ECFP4 ~0.7; FCFP4 ~0.6.

## Reconciliation: Cliffs Across Methods

If a pair flags as an activity cliff in ECFP4 but not in MCS-based, the local environment differs at a single key atom (e.g., -F to -OH). If it flags in MCS but not ECFP4, a remote substituent affects activity (e.g., para vs meta). Use both views to localize the SAR.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `BulkTanimotoSimilarity` returns ints | Bit-vector with single bit per atom | Already correct; ints are bit counts |
| Tanimoto > 1 | Wrong coefficient (Tversky alpha+beta != 1, count vector) | Use bit vector + Tanimoto |
| Cluster centroids change between runs | Order-dependent Butina | Sort molecules first or use stable random seed for tie-breaking |
| MaxMinPicker returns first N inputs | All-zero initial similarity matrix | Seed picker explicitly: `picker.LazyBitVectorPick(fps, n_lib, n_pick, seed=42)` |
| Activity cliff "false positives" | Bit-collisions inflate similarity | Use sparse Morgan or compare canonical SMILES for exact ID |
| Diverse subset has duplicates | Standardization not applied | Canonicalize via `chemoinformatics/molecular-standardization` first |
| Tanimoto incompatible with neural fingerprint | Continuous-valued fingerprint | Use cosine or sklearn `cdist` with `'cosine'` metric |

## References

- Bajorath, *Nat. Rev. Drug Discov.* 1:882 (2002) -- molecular similarity principle review.
- Maggiora et al., *J. Med. Chem.* 57:3186 (2014) -- molecular similarity in drug discovery.
- Stumpfe & Bajorath, *J. Med. Chem.* 55:2932 (2012) -- activity cliffs.
- Tversky, *Psychol. Rev.* 84:327 (1977) -- features of similarity (Tversky coefficient).
- Probst & Reymond, *J. Cheminformatics* 10:66 (2018) -- MHFP6 MinHash fingerprint.
- O'Boyle et al., *J. Cheminformatics* 8:36 (2016) -- comparing fingerprints for similarity.
- Butina, *J. Chem. Inf. Comput. Sci.* 39:747 (1999) -- Taylor-Butina clustering.

## Related Skills

- chemoinformatics/molecular-descriptors - Generate fingerprints for similarity
- chemoinformatics/molecular-standardization - Canonicalize before comparing
- chemoinformatics/substructure-search - SMARTS pattern-based searching
- chemoinformatics/scaffold-analysis - Scaffold-based similarity (Bemis-Murcko, MMPA)
- chemoinformatics/shape-similarity - 3D shape similarity (USRCAT, ROCS)
- machine-learning/biomarker-discovery - ML on similarity features
