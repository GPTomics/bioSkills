---
name: bio-molecular-standardization
description: Standardizes molecular structures using the ChEMBL structure pipeline for normalization and parent selection plus RDKit rdMolStandardize for explicit custom steps such as tautomer canonicalization, salt/solvent stripping, charge handling, stereochemistry handling, mixture selection, and isotope normalization. Explicitly compares ChEMBL, canSARchem, RDKit, and PubChem standardization choices. Use when preparing libraries for QSAR training, joining datasets across sources, deduplicating compound collections, or building canonical compound registries.
tool_type: python
primary_tool: RDKit
---

## Version Compatibility

Reference examples tested with: RDKit 2024.09+ and chembl_structure_pipeline 1.2+. MolVS 0.1.1 is a legacy package; use RDKit's maintained `rdMolStandardize` module for custom pipelines.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Molecular Standardization

Convert raw molecular structures into a consistent form for ML training data, deduplication, registry, and cross-database joining. Skipping standardization can create data leakage when alternate representations of one compound enter different splits, distort QSAR inputs, and cause database join misses. The ChEMBL structure pipeline (Bento et al. 2020) is built on RDKit and applies ChEMBL-specific normalization and parent-selection rules. canSARchem (Dolciami et al. 2022) adds canonical-tautomer selection before parent extraction. RDKit's maintained `rdMolStandardize` module provides primitives for building an explicit custom pipeline.

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
| 7. Isotope normalization | Explicitly set selected atom isotope labels to 0 | Remove 13C, 2H labels | Tracer studies; preserve labels when scientifically meaningful |
| 8. Output canonicalization | `Chem.MolToSmiles(canonical=True)` | Canonical SMILES + InChIKey | Round-trip stability |

## Pipeline Reconciliation

| Pipeline | Origin | Tautomer canonicalization | Salt definition | Use case |
|----------|--------|---------------------------|-----------------|----------|
| ChEMBL pipeline | EBI ChEMBL | Not performed by `standardize_mol` or `get_parent_mol` | ChEMBL salt list (extensive) | ChEMBL-compatible registration |
| canSARchem | ICR Cancer Research UK | Canonical tautomer BEFORE parent extraction | Extended salt list | Cancer drug discovery |
| PubChem (OpenEye) | NIH NCBI | OpenEye QUACPAC tautomer | PubChem salt list | Bioassay data, large-scale |
| RDKit rdMolStandardize default | Greg Landrum | RDKit TautomerEnumerator | RDKit default | General purpose, open source |

**Key difference (canSARchem vs ChEMBL):**
- ChEMBL standardizes the representation and extracts a parent, but does not canonicalize tautomers.
- canSARchem canonicalizes the tautomer before parent extraction.

This difference matters when alternate tautomeric inputs must be registered as one parent. Do not describe ChEMBL output as tautomer-canonical unless an explicit tautomer step is added and documented.

## ChEMBL Structure Pipeline (Reference Implementation)

ChEMBL's standardization is the most widely-used reference. The Python package `chembl_structure_pipeline` exposes the validated pipeline.

**Goal:** Apply the industry-reference ChEMBL standardization pipeline to a SMILES.

**Approach:** Parse SMILES with RDKit, run `standardize_mol` (sanitize, normalize, and standardize charges), then `get_parent_mol` (strip salts/counter-ions), and emit canonical SMILES. Add `rdMolStandardize.TautomerEnumerator` separately only when the project requires tautomer canonicalization.

```python
from chembl_structure_pipeline import standardize_mol, get_parent_mol
from rdkit import Chem

def chembl_pipeline(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None, 'parse_failure'
    standardized = standardize_mol(mol)
    parent, exclude = get_parent_mol(standardized)
    if exclude:
        return None, 'excluded_by_chembl'
    return Chem.MolToSmiles(parent), 'ok'
```

**`standardize_mol`:** sanitize, normalize functional groups, and standardize charges; returns one RDKit molecule.

**`get_parent_mol`:** strip salts/counter-ions and choose the parent; returns `(parent_mol, exclude_flag)`.

Output: canonical SMILES of the selected parent after the ChEMBL transformations, or an explicit `excluded_by_chembl` status when the parent carries ChEMBL's exclusion flag. Neutralizable acid/base sites may be normalized, but permanent or otherwise non-removable charges can remain; do not assume every emitted parent is neutral.

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

    uncharger = rdMolStandardize.Uncharger(canonicalOrder=True)
    mol = uncharger.uncharge(mol)

    enumerator = rdMolStandardize.TautomerEnumerator()
    mol = enumerator.Canonicalize(mol)

    if not keep_isotopes:
        for atom in mol.GetAtoms():
            atom.SetIsotope(0)

    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return Chem.MolToSmiles(mol)
```

**`canonicalOrder=True`** makes the uncharger choose neutralization sites in canonical order when more than one equivalent site is available. It does not itself decide whether a permanent charge is retained; inspect charge-sensitive structures and keep `force=False` unless a documented policy requires otherwise.

## Salt Stripping Edge Cases

| Salt form | Action | Example |
|-----------|--------|---------|
| Mono-salt | Strip counter-ion | `[Na+].CC(=O)[O-]` -> `CC(=O)O` |
| Di-salt | Strip both | `[Na+].[Na+].CC(=O)[O-].CC(=O)[O-]` -> `CC(=O)O` |
| Mixed salt | Largest organic fragment | `CCO.CC(=O)O` -> `CCO` (or CC(=O)O depending on rule) |
| Co-crystal | Hardest case | `CC(=O)O.CCOC(C)=O` -- both organic; default returns largest |
| Hydrate | Strip waters | `CC(=O)O.O` -> `CC(=O)O` |
| Solvate | Strip solvents | `CC(=O)O.CO` -> `CC(=O)O` |
| Quaternary ammonium | Preserve charge | `[N+](C)(C)(C)C` (permanent charge; do NOT neutralize) |

**`LargestFragmentChooser(preferOrganic=True)`** prefers organic fragments over inorganic counter-ions even if smaller; for co-crystals, default rule picks largest organic fragment.

## Tautomer Canonicalization (debated)

Tautomer canonicalization is the most controversial standardization step. There is no universally-correct canonical tautomer for many drug-like molecules.

| Tautomer pair | Why the policy matters |
|---------------|------------------------|
| Keto/enol | Canonicalization can select a representation different from the experimentally relevant bound or solution form |
| Lactam/lactim | Heterocycle scoring rules and toolkit versions may choose different representatives |
| Amidine/iminol | Proton placement changes donor/acceptor annotations and downstream matching |
| Phenol/keto (e.g., naphthol/naphthalenone) | Aromaticity and functional-group perception can change with the selected representation |
| 2H-pyrazole / 1H-pyrazole | Nitrogen identity and donor/acceptor assignments depend on proton placement |

Treat the enumerator's canonical result as a reproducible representation chosen by its configured scoring rules, not as a prediction of the dominant tautomer in vivo. Record the RDKit version and any custom transforms or scoring changes.

**Practical rules:**
- Always apply consistent canonicalization across train + test for ML
- For prospective prediction, predict for both tautomers if disagreement could matter
- For library deduplication, canonical tautomer is the standard answer
- For docking, use an ionization-aware preparation workflow. For Open Babel, the documented CLI is `obabel input.sdf -O output.sdf -p 7.4`; validate generated states because its rule-based protonation is not a substitute for project-specific pKa analysis.

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
- Explicit stereo with `@` / `\` / `/` -> preserved
- Wedge bonds in SDF -> re-perceived from 3D coords if present
- Ambiguous stereo (no markers) -> left as-is, marked as undefined
- Racemic (explicit "rac") -> keep as racemate

For ML, remove stereochemistry only when the endpoint, data curation, and model representation justify treating stereoisomers as equivalent; record that policy and test its effect. For docking and FEP, preserve the intended stereoisomer and reject unintended stereo changes.

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
            mol = standardize_mol(mol)
            mol, exclude = get_parent_mol(mol)
            if exclude:
                continue
            standardized.append({
                'smiles': Chem.MolToSmiles(mol),
                'inchikey': Chem.MolToInchiKey(mol),
                'activity': row[activity_col],
            })
        except Exception:
            continue

    df_std = pd.DataFrame(standardized)
    if df_std.empty:
        return pd.DataFrame(columns=['inchikey', 'smiles', 'activity', 'n_replicates'])
    df_std = df_std.groupby('inchikey').agg(
        smiles=('smiles', 'first'),
        activity=('activity', 'mean'),
        n_replicates=('activity', 'count'),
    ).reset_index()
    return df_std
```

Standard InChIKey may collapse some mobile-hydrogen tautomer representations, but this is not a substitute for an explicitly chosen tautomer policy. Replicate count signals measurement reliability.

## Per-Tool Failure Modes

### ChEMBL pipeline -- inorganic salt fails

**Trigger:** Molecule is genuinely an inorganic salt (e.g., NaCl, K2SO4).

**Mechanism:** `get_parent_mol` chooses largest organic; falls back to largest fragment for fully inorganic.

**Symptom:** Returns the salt itself (not a drug).

**Fix:** Pre-filter to compounds with ≥1 carbon atom.

### Uncharger -- charge-state policy mismatch

**Trigger:** A molecule combines a non-removable charge, such as quaternary ammonium, with other neutralizable sites, or the desired physiological ionization state differs from a structure-normalization rule.

**Mechanism:** `Uncharger` adds or removes hydrogens from neutralizable acids and bases. It cannot remove a permanent charge that has no corresponding hydrogen edit; by default it may preserve an opposite neutralizable charge when a non-removable charge is present so that the total charge remains balanced. `force=True` instead neutralizes all sites that can be neutralized even if the remaining permanent charge leaves a nonzero total charge.

**Symptom:** The permanent charge remains, but other sites or the total charge differ from the protonation state intended for docking or modeling.

**Fix:** Choose `force` according to the documented total-charge policy, keep `force=False` when balanced countercharges should be preserved, and inspect/prepare physiological protonation states separately.

### Tautomer enumerator -- combinatorial explosion

**Trigger:** Molecule with many tautomerizable groups (polyhydroxylated heterocycle).

**Mechanism:** `TautomerEnumerator.Enumerate` generates all possible tautomers; can produce thousands.

**Symptom:** OOM or hour-long compute on single molecule.

**Fix:** Use `Canonicalize` when only the configured canonical representation is needed. Before `Enumerate`, call `enumerator.SetMaxTransforms(limit)` (and, when appropriate, `SetMaxTautomers(limit)`) to cap the search.

### Legacy MolVS -- import or compatibility failure

**Trigger:** Code still using legacy `from molvs import Standardizer`.

**Mechanism:** The standalone MolVS package is legacy and may not support current Python/RDKit versions. RDKit's maintained `rdMolStandardize` module remains available.

**Symptom:** ImportError or AttributeError on newer RDKit.

**Fix:** Migrate deliberately to `from rdkit.Chem.MolStandardize import rdMolStandardize`; compare outputs because RDKit functions are not drop-in aliases for every MolVS workflow.

### Round-trip InChIKey mismatch

**Trigger:** Records were processed with different standardization settings or entered in different salt, charge, isotope, stereo, or tautomer forms.

**Mechanism:** The pipelines did not apply the same explicitly versioned transformations before identity generation.

**Symptom:** Apparently equivalent records produce different InChIKeys, or an expected database join fails.

**Fix:** Record and apply the same toolkit version, standardization stages, tautomer policy, and InChI options to both datasets; compare full standardized structures when results still differ.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| ImportError from standalone `molvs` | Legacy package incompatible with current environment | Use maintained `rdkit.Chem.MolStandardize.rdMolStandardize` APIs and validate output |
| `standardize_mol` raises or input parsing returns `None` | Invalid or unsanitizable input | Capture the exception/input index and inspect sanitization deliberately; do not silently accept a partially sanitized structure |
| Stripped wrong fragment | LargestFragmentChooser ambiguity | Manually inspect; consider custom logic |
| Tautomer differs between datasets | Different tautomer rules or toolkit versions | Pin and record the same `TautomerEnumerator` settings and version |
| Unexpected charge distribution with permanent ions | `Uncharger` total-charge policy does not match the intended protonation workflow | Review non-removable and neutralizable sites; choose `force` deliberately and prepare physiological states separately |
| Same InChIKey for apparently different records | Standard-InChI normalization or a rare hash collision | Compare full InChI and standardized structures; InChIKey has no longer form |
| Pipeline slow on large library | Per-molecule Python overhead | Process independent molecules in validated chunks or worker processes; `chembl_structure_pipeline` itself is a per-molecule API |

## References

- Bento et al., *J. Cheminformatics* 12:51 (2020) -- ChEMBL structure pipeline. https://doi.org/10.1186/s13321-020-00456-1
- Dolciami et al., *J. Cheminformatics* 14:28 (2022) -- canSARchem registration pipeline. https://doi.org/10.1186/s13321-022-00606-7
- Hähnke et al., *J. Cheminformatics* 10:36 (2018) -- PubChem standardization. https://doi.org/10.1186/s13321-018-0293-8
- RDKit, `rdkit.Chem.MolStandardize.rdMolStandardize` API documentation. https://www.rdkit.org/docs/source/rdkit.Chem.MolStandardize.rdMolStandardize.html
- Open Babel, official `obabel` documentation -- pH-dependent hydrogen-addition CLI. https://openbabel.org/docs/Command-line_tools/babel.html

## Related Skills

- chemoinformatics/molecular-io - Parse molecules before standardizing
- chemoinformatics/molecular-descriptors - Apply descriptors to standardized molecules
- chemoinformatics/similarity-searching - Standardize before comparing
- chemoinformatics/substructure-search - Standardize before SMARTS matching
- chemoinformatics/qsar-modeling - Mandatory upstream for QSAR
