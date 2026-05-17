---
name: bio-molecular-standardization
description: Standardizes molecular structures using ChEMBL chembl_structure_pipeline and RDKit rdMolStandardize covering sanitization, salt/solvent stripping, neutralization, tautomer canonicalization, stereochemistry standardization, mixture handling, and isotope normalization. Explicitly compares ChEMBL pipeline, canSARchem, and PubChem standardization choices. Use when preparing libraries for QSAR training, joining datasets across sources, deduplicating compound collections, or building canonical compound registries.
tool_type: python
primary_tool: RDKit
---

## Version Compatibility

Reference examples tested with: RDKit 2024.09+, chembl_structure_pipeline 1.2+, MolVS 0.1.1 (legacy reference only -- rdMolStandardize is current).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Molecular Standardization

Convert raw molecular structures into a single canonical form for ML training data, deduplication, registry, and cross-database joining. Standardization is the single most underrated upstream step: skipping it causes silent ML data leakage (training and test compounds with different tautomers count as separate), bogus QSAR predictions, and database join misses. The ChEMBL pipeline (Bento 2020) and canSARchem (Ravi 2022) are the two industry references; canSARchem extends ChEMBL with canonical-tautomer-before-parent extraction. RDKit's `rdMolStandardize` implements ChEMBL-equivalent logic in C++ (the older `MolVS` Python implementation was deprecated Q1 2024).

For format-level I/O and aromaticity perception, see `chemoinformatics/molecular-io`. For descriptor calculation after standardization, see `chemoinformatics/molecular-descriptors`.

## Standardization Pipeline Stages

| Stage | RDKit Tool | Operation | Common errors caught |
|-------|-----------|-----------|----------------------|
| 1. Sanitization | `Chem.SanitizeMol` | Kekulize, assign aromaticity, fix valences | Wrong valence on N/O |
| 2. Salt stripping | `rdMolStandardize.FragmentRemover` or `LargestFragmentChooser` | Remove counterions | Cl-, Na+, K+, OH- |
| 3. Mixture choice | `LargestFragmentChooser` | Pick parent fragment | Co-crystals, hydrates |
| 4. Charge neutralization | `Uncharger` | Neutralize while preserving net charge | Permanent charges preserved (quaternary N+) |
| 5. Tautomer canonicalization | `TautomerEnumerator.Canonicalize` | Pick canonical tautomer | Keto/enol; amide/imidate |
| 6. Stereo standardization | `Chem.AssignStereochemistry` | Consistent stereo descriptors | Lost wedges, ambiguous R/S |
| 7. Isotope normalization | manual or `MolToSmiles(isomericSmiles=False)` | Remove 13C, 2H labels | Tracer studies |
| 8. Output canonicalization | `Chem.MolToSmiles(canonical=True)` | Canonical SMILES + InChIKey | Round-trip stability |

## Pipeline Reconciliation

| Pipeline | Origin | Tautomer canonicalization | Salt definition | Use case |
|----------|--------|---------------------------|-----------------|----------|
| ChEMBL pipeline | EBI ChEMBL | Pre-rdMolStandardize legacy; now uses rdMolStandardize | ChEMBL salt list (extensive) | Drug-like compounds, FDA approvals |
| canSARchem | ICR Cancer Research UK | Canonical tautomer BEFORE parent extraction | Extended salt list | Cancer drug discovery |
| PubChem (OpenEye) | NIH NCBI | OpenEye QUACPAC tautomer | PubChem salt list | Bioassay data, large-scale |
| RDKit rdMolStandardize default | Greg Landrum | RDKit TautomerEnumerator | RDKit default | General purpose, open source |

**Key difference (canSARchem vs ChEMBL):**
- ChEMBL: parent first, then canonical tautomer of parent
- canSARchem: canonical tautomer first, then parent

For 95% of drug-like molecules these produce identical results. For tautomer-ambiguous molecules (amide/imidate, ketoenol, lactam/lactim), the order matters; canSARchem produces more stable canonical forms.

## ChEMBL Structure Pipeline (Reference Implementation)

ChEMBL's standardization is the most widely-used reference. The Python package `chembl_structure_pipeline` exposes the validated pipeline.

**Goal:** Apply the industry-reference ChEMBL standardization pipeline to a SMILES.

**Approach:** Parse SMILES with RDKit, run `standardize_mol` (sanitize + uncharge + normalize + canonical tautomer), then `get_parent_mol` (strip salts/counter-ions), and emit canonical SMILES.

```python
from chembl_structure_pipeline import standardize_mol, get_parent_mol
from rdkit import Chem

def chembl_pipeline(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None, 'parse_failure'
    standardized, _ = standardize_mol(mol)
    parent, _ = get_parent_mol(standardized)
    return Chem.MolToSmiles(parent), 'ok'
```

**`standardize_mol`:** sanitize + uncharge + normalize functional groups + canonicalize tautomers.

**`get_parent_mol`:** strip salts/counter-ions; choose largest fragment.

Output: canonical SMILES of the parent (free acid/free base, neutral form).

## Full Standardization with rdMolStandardize

For more granular control or non-ChEMBL workflows.

**Goal:** Execute each standardization step explicitly to control salt stripping, charge handling, tautomer canonicalization, and isotope normalization.

**Approach:** Run the 8-stage pipeline (sanitize, largest fragment, normalize, uncharge, tautomer canonicalize, isotope strip, stereo standardize, canonical SMILES) sequentially with `rdMolStandardize` primitives.

```python
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

def full_standardize(smi, keep_isotopes=False):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None

    Chem.SanitizeMol(mol)

    largest = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)
    mol = largest.choose(mol)

    normalizer = rdMolStandardize.Normalizer()
    mol = normalizer.normalize(mol)

    uncharger = rdMolStandardize.Uncharger(canonicalOrdering=True)
    mol = uncharger.uncharge(mol)

    enumerator = rdMolStandardize.TautomerEnumerator()
    mol = enumerator.Canonicalize(mol)

    if not keep_isotopes:
        for atom in mol.GetAtoms():
            atom.SetIsotope(0)

    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return Chem.MolToSmiles(mol)
```

**`canonicalOrdering=True`** ensures the uncharger produces the same result regardless of atom ordering in input -- critical for stable canonical output.

## Salt Stripping Edge Cases

| Salt form | Action | Example |
|-----------|--------|---------|
| Mono-salt | Strip counter-ion | `[Na+].CC(=O)[O-]` → `CC(=O)O` |
| Di-salt | Strip both | `[Na+].[Na+].CC(=O)[O-].CC(=O)[O-]` → `CC(=O)O` |
| Mixed salt | Largest organic fragment | `CCO.CC(=O)O` → `CCO` (or CC(=O)O depending on rule) |
| Co-crystal | Hardest case | `CC(=O)O.CCOC(C)=O` -- both organic; default returns largest |
| Hydrate | Strip waters | `CC(=O)O.O` → `CC(=O)O` |
| Solvate | Strip solvents | `CC(=O)O.CO` → `CC(=O)O` |
| Quaternary ammonium | Preserve charge | `[N+](C)(C)(C)C` (permanent charge; do NOT neutralize) |

**`LargestFragmentChooser(preferOrganic=True)`** prefers organic fragments over inorganic counter-ions even if smaller; for co-crystals, default rule picks largest organic fragment.

## Tautomer Canonicalization (debated)

Tautomer canonicalization is the most controversial standardization step. There is no universally-correct canonical tautomer for many drug-like molecules.

| Tautomer pair | Default canonical | Issue |
|---------------|-------------------|-------|
| Keto/enol | Keto preferred | Most kinase ATP-mimetic enols destabilize on canonicalization |
| Lactam/lactim | Lactam preferred | Some natural products (rifampin) are inherently lactim |
| Amidine/iminol | Amidine preferred | Some bioactive amidines convert |
| Phenol/keto (e.g., naphthol/naphthalenone) | Phenol preferred | Some quinone-form pharmaceuticals reverted |
| 2H-pyrazole / 1H-pyrazole | 1H-pyrazole | Both equally stable in vivo |

**Practical rules:**
- Always apply consistent canonicalization across train + test for ML
- For prospective prediction, predict for both tautomers if disagreement could matter
- For library deduplication, canonical tautomer is the standard answer
- For docking, use ionization-aware preparation (`epik` from Schrödinger or `Open Babel pkBABEL`)

```python
from rdkit.Chem.MolStandardize import rdMolStandardize

def canonical_tautomer(smi):
    mol = Chem.MolFromSmiles(smi)
    enumerator = rdMolStandardize.TautomerEnumerator()
    canon = enumerator.Canonicalize(mol)
    return Chem.MolToSmiles(canon)
```

## Stereochemistry Standardization

```python
from rdkit import Chem

def standardize_stereo(mol, remove_undefined=False):
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    if remove_undefined:
        Chem.RemoveStereochemistry(mol)
    return mol
```

**Cases:**
- Explicit stereo with `@` / `\` / `/` → preserved
- Wedge bonds in SDF → re-perceived from 3D coords if present
- Ambiguous stereo (no markers) → left as-is, marked as undefined
- Racemic (explicit "rac") → keep as racemate

For ML, often drop stereo entirely (`Chem.RemoveStereochemistry(mol)`) since most QSAR endpoints are not stereo-specific. For docking and FEP, preserve stereo always.

## Standardization for ML Training (avoiding data leakage)

**Goal:** Build a standardized + deduplicated training set with replicate-averaged activity for QSAR or ADMET model training.

**Approach:** Standardize every SMILES through the ChEMBL pipeline, compute InChIKey as canonical identity, group by InChIKey, and mean-aggregate activities; report replicate count for confidence weighting.

```python
import pandas as pd
from chembl_structure_pipeline import standardize_mol, get_parent_mol

def prepare_qsar_data(df, smiles_col='smiles', activity_col='pIC50'):
    standardized = []
    for i, row in df.iterrows():
        mol = Chem.MolFromSmiles(row[smiles_col])
        if mol is None:
            continue
        try:
            mol, _ = standardize_mol(mol)
            mol, _ = get_parent_mol(mol)
            standardized.append({
                'smiles': Chem.MolToSmiles(mol),
                'inchikey': Chem.MolToInchiKey(mol),
                'activity': row[activity_col],
            })
        except Exception:
            continue

    df_std = pd.DataFrame(standardized)
    df_std = df_std.groupby('inchikey').agg(
        smiles=('smiles', 'first'),
        activity=('activity', 'mean'),
        n_replicates=('activity', 'count'),
    ).reset_index()
    return df_std
```

Deduplication by InChIKey collapses tautomer-equivalent compounds. Replicate count signals measurement reliability.

## Per-Tool Failure Modes

### ChEMBL pipeline -- inorganic salt fails

**Trigger:** Molecule is genuinely an inorganic salt (e.g., NaCl, K2SO4).

**Mechanism:** `get_parent_mol` chooses largest organic; falls back to largest fragment for fully inorganic.

**Symptom:** Returns the salt itself (not a drug).

**Fix:** Pre-filter to compounds with ≥1 carbon atom.

### Uncharger -- removes critical charge

**Trigger:** Quaternary ammonium (permanent positive) or sulfonate at physiological pH.

**Mechanism:** Default uncharger attempts to neutralize without distinguishing permanent vs pH-dependent charges.

**Symptom:** Permanently charged ligands neutralized; structure incorrect for downstream docking.

**Fix:** Use `Uncharger(canonicalOrdering=True, force=False)`; manually inspect borderline cases.

### Tautomer enumerator -- combinatorial explosion

**Trigger:** Molecule with many tautomerizable groups (polyhydroxylated heterocycle).

**Mechanism:** `TautomerEnumerator.Enumerate` generates all possible tautomers; can produce thousands.

**Symptom:** OOM or hour-long compute on single molecule.

**Fix:** Use `Canonicalize` (returns single canonical) instead of `Enumerate`; for `Enumerate`, cap `maxTransforms` parameter.

### MolVS deprecated -- ImportError on Python 3.12+

**Trigger:** Code still using legacy `from molvs import Standardizer`.

**Mechanism:** RDKit MolStandardize Python implementation removed Q1 2024.

**Symptom:** ImportError or AttributeError on newer RDKit.

**Fix:** Migrate to `from rdkit.Chem.MolStandardize import rdMolStandardize`; methods renamed (e.g., `standardize` → `Cleanup`).

### Round-trip InChIKey mismatch

**Trigger:** Compound canonicalized to different tautomer per run.

**Mechanism:** RDKit tautomer canonicalization depends on atom ordering for very symmetric molecules.

**Symptom:** Re-running standardization yields different InChIKey.

**Fix:** Set `canonicalOrdering=True` in Uncharger; sort atoms via canonical SMILES first.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| ImportError on `MolStandardize` | Python MolStandardize deprecated | Use `from rdkit.Chem.MolStandardize import rdMolStandardize` |
| `standardize_mol` returns None | Sanitize failure on input | Try `Chem.MolFromSmiles(smi, sanitize=False)` first |
| Stripped wrong fragment | LargestFragmentChooser ambiguity | Manually inspect; consider custom logic |
| Tautomer differs between runs | Atom-order-dependent | Set `canonicalOrdering=True`; sort atoms |
| Charge lost on quaternary N | Aggressive neutralization | Use `force=False` |
| InChIKey collisions across "different" mols | Same canonical InChI but different stereo / tautomer | Use longer InChIKey (or full InChI) |
| Pipeline slow on large library | Per-mol Python overhead | Use `chembl_structure_pipeline` (vectorized) or process in chunks |

## References

- Bento et al., *J. Cheminformatics* 12:51 (2020) -- ChEMBL structure pipeline.
- Ravi et al., *J. Cheminformatics* 14:36 (2022) -- canSARchem registration pipeline.
- Hähnke et al., *J. Cheminformatics* 10:36 (2018) -- PubChem standardization.
- Landrum, RDKit Blog (2020) -- new tautomer canonicalization code.

## Related Skills

- chemoinformatics/molecular-io - Parse molecules before standardizing
- chemoinformatics/molecular-descriptors - Apply descriptors to standardized molecules
- chemoinformatics/similarity-searching - Standardize before comparing
- chemoinformatics/substructure-search - Standardize before SMARTS matching
- chemoinformatics/qsar-modeling - Mandatory upstream for QSAR
