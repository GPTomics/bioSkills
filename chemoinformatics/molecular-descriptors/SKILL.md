---
name: bio-molecular-descriptors
description: Calculates molecular fingerprints (ECFP/Morgan, FCFP, MACCS, RDKit, AtomPair, TopologicalTorsion, Avalon, MAP4, MHFP6) and physicochemical descriptors (Lipinski, QED, TPSA, Crippen LogP, 3D shape) with explicit choice tables, bit vs count semantics, and partial-charge model selection. Use when featurizing molecules for similarity, QSAR, virtual screening, or ML, or selecting the correct fingerprint for a chemotype-aware task.
tool_type: python
primary_tool: RDKit
---

## Version Compatibility

Reference examples tested with: RDKit 2024.09+, numpy 1.26+, pandas 2.2+, mapchiral 0.1+ (MAP4), mhfp 1.9+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Molecular Descriptors

Featurize molecules for similarity search, QSAR, virtual screening, or ML. The fingerprint or descriptor choice is **chemotype-aware**: ECFP4 dominates drug-like organic similarity, AtomPair and TopologicalTorsion outperform for scaffold hopping, MAP4/MHFP6 win on metabolomics-scale chemical diversity, and 3D conformer-based descriptors are essential when shape and stereochemistry matter.

For canonicalization before featurization, see `chemoinformatics/molecular-standardization`. For 3D-only descriptors, see `chemoinformatics/conformer-generation`.

## Fingerprint Taxonomy

| Fingerprint | Type | Radius/Path | Bits | Use case | Fails when |
|-------------|------|-------------|------|----------|------------|
| Morgan (ECFP) | Circular | r=2 (ECFP4), r=3 (ECFP6) | 2048 typical | Drug-like similarity, ML default | Loses long-range topology; bit collisions at low nBits |
| FCFP | Functional Morgan | r=2 default | 2048 | Pharmacophore-aware similarity | Same caveats as ECFP; less specific |
| MACCS | Substructure key | 166 fixed bits | 167 | Quick fingerprint, drug-likeness | Too sparse for large diverse libraries |
| RDKit FP | Path-based | linear paths up to 7 atoms | 2048 | RDKit-native ECFP alternative | Drug-like only; not optimal for scaffold hopping |
| AtomPair | Pair + topological distance | All atom pairs | 2048 | Scaffold hopping; flexible mol | Slower than ECFP; harder to interpret |
| TopologicalTorsion | 4-atom torsion | All TT | 2048 | Scaffold hopping; less hit-rate | Like AP, slower than ECFP |
| Avalon | Substructure + atom pairs | Mixed | 512/1024 | Fast similarity | Less standard; older |
| MAP4 (MinHashed atom-pair) | MinHash atom-pair | r=1,2 | 1024/2048 | Biological + metabolite diversity | Library required (mapchiral); slower hash |
| MHFP6 (MinHash) | MinHash ECFP-like | r=3 (diam 6) | 2048 | Big-data nearest-neighbor (Annoy) | Different distance (Jaccard on MinHash) |
| Pharm2D | 2D pharmacophore | feature pairs/triplets | sparse | Pharmacophore search | Sparse, slower |

**Decision:** For drug-like similarity ranking, use **ECFP4 2048 bit**; established baseline, fast, well-understood. For diverse libraries (>1M compounds, metabolomics, peptides), **MHFP6** outperforms ECFP4 on analog recovery (Probst & Reymond 2018). For scaffold hopping, **AtomPair** beats ECFP4 on retrospective benchmarks but loses on retrospective single-target.

## Bit vs Count Vectors

| Form | Use | Library impact |
|------|-----|----------------|
| Bit (0/1) | Tanimoto similarity, BulkTanimotoSimilarity, RDKit fingerprint folding | Standard for similarity |
| Count (integer) | Some ML methods, RF on counts, neural fingerprints | Loses bit-level fast operations; richer signal |
| Sparse (dict) | Direct chemical interpretation (which fragments at which atoms) | Use for SHAP / atomic attribution |

```python
from rdkit import Chem
from rdkit.Chem import AllChem

mol = Chem.MolFromSmiles('CCO')

ecfp4_bit = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
ecfp4_count = AllChem.GetHashedMorganFingerprint(mol, radius=2, nBits=2048)
ecfp4_sparse = AllChem.GetMorganFingerprint(mol, radius=2)
```

## Morgan / ECFP Radius Math

ECFP-X notation: X is the **diameter** in bonds. RDKit's `radius` parameter is half of X.

| Notation | RDKit radius | Diameter | Captures |
|----------|--------------|----------|----------|
| ECFP0 | 0 | 0 | Atom identity only |
| ECFP2 | 1 | 2 | Atom + immediate neighbors |
| ECFP4 | 2 | 4 | Atom + 2-bond environment |
| ECFP6 | 3 | 6 | Atom + 3-bond environment |

**Trade-off:** Larger radius captures more specific local environment but inflates bit-collision rate at fixed nBits. For QSAR with <10k compounds, **ECFP4 2048** is the established default (Rogers & Hahn 2010; MoleculeNet benchmarks Wu 2018). For large libraries (>1M), use **nBits=4096** or unhashed sparse representation to reduce ~1-5% bit-collision rate (O'Boyle 2016).

## FCFP vs ECFP

FCFP (Functional-Class) uses atom invariants based on pharmacophore role (donor, acceptor, hydrophobe, aromatic, halogen, basic, acidic) instead of atom identity. Trades atom-specificity for functional-equivalence.

```python
ecfp4 = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, useFeatures=False)
fcfp4 = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, useFeatures=True)
```

**When to use FCFP4:** Scaffold-hopping campaigns, pharmacophore-driven similarity, cross-target activity prediction.

**When to use ECFP4:** Within-series QSAR, lead optimization, when chemotype identity matters.

## 3D Descriptors and Conformer Dependence

Conformer-dependent descriptors (asphericity, eccentricity, principal moments of inertia, RDF) require a generated 3D structure. A **single conformer** is rarely sufficient: descriptor variance across the conformer ensemble can exceed the descriptor signal.

**Goal:** Compute 3D shape descriptors over a conformer ensemble rather than from a single (possibly unrepresentative) conformer.

**Approach:** Add explicit hydrogens, embed N conformers with ETKDGv3, MMFF-optimize them all, then evaluate the descriptor across each conformer for downstream averaging.

```python
from rdkit.Chem import AllChem, Descriptors3D

mol = Chem.MolFromSmiles('CCCCO')
mol = Chem.AddHs(mol)

params = AllChem.ETKDGv3()
params.randomSeed = 42
n = AllChem.EmbedMultipleConfs(mol, numConfs=20, params=params)
AllChem.MMFFOptimizeMoleculeConfs(mol)

asphericities = [Descriptors3D.Asphericity(mol, confId=c) for c in range(n)]
```

**Decision:** For QSAR / ML, compute over a conformer ensemble (n=20-100) and report mean or Boltzmann-weighted average. Single-conformer 3D descriptors are unreliable.

## Partial Charge Methods

| Method | Software | Cost | Accuracy | Use for |
|--------|----------|------|----------|---------|
| Gasteiger-Marsili | RDKit, Open Babel | <0.1s/mol | Empirical, rough | AutoDock Vina, fast screening |
| MMFF94 | RDKit | 0.1s/mol | Force-field consistent | MMFF energy, conformer ranking |
| AM1-BCC | antechamber (AmberTools) | ~10s/mol | Semi-empirical | MD setup, FEP, GAFF |
| RESP | psi4, Gaussian | minutes/mol | DFT ESP-fitted | High-accuracy MD, FEP |
| OpenFF Recharge | openff-recharge | seconds | DFT-derived but cached | OpenFF / SAGE setup |
| ABCG2 | Open Babel | <1s | Improved empirical | Modern Vina, AutoDock-GPU |

```python
from rdkit.Chem import AllChem

AllChem.ComputeGasteigerCharges(mol)
for atom in mol.GetAtoms():
    print(atom.GetIdx(), atom.GetPropsAsDict().get('_GasteigerCharge', None))
```

**Critical:** Charge method must match downstream. Gasteiger charges in an AMBER MD run violate the assumptions of the protein force field.

## MAP4 and MHFP6 for Diverse Libraries

For libraries spanning drug-like + natural products + peptides + metabolites, ECFP4 saturates Tanimoto similarity (most pairs report 0.1-0.3, hard to rank). MAP4 and MHFP6 use MinHash + atom-pair / circular substructures and discriminate better.

```python
from mhfp.encoder import MHFPEncoder

encoder = MHFPEncoder(2048)
mhfp6 = encoder.encode(mol, radius=3)
```

MHFP6 distance is Jaccard on MinHash, not standard Tanimoto. Use `MHFPEncoder.distance(fp1, fp2)`.

## Physicochemical Descriptors

| Descriptor | Source | Range | Drug-like cutoff |
|------------|--------|-------|-------------------|
| MolWt | RDKit `Descriptors.MolWt` | ~50-2000 Da | <=500 (Lipinski) |
| MolLogP (Crippen) | RDKit `Descriptors.MolLogP` | -5 to 8 | <=5 (Lipinski) |
| HBD | `Lipinski.NumHDonors` | 0-10 | <=5 (Lipinski) |
| HBA | `Lipinski.NumHAcceptors` | 0-15 | <=10 (Lipinski) |
| TPSA | `Descriptors.TPSA` (Ertl) | 0-200 A^2 | <=140 (Veber oral); <=90 (BBB+) |
| RotBonds | `Lipinski.NumRotatableBonds` | 0-15 | <=10 (Veber) |
| AromaticRings | `Lipinski.NumAromaticRings` | 0-6 | <=3-4 (Ritchie-Macdonald aromatic ring count) |
| HeavyAtoms | `Descriptors.HeavyAtomCount` | <=50 (lead-like) | |
| FractionCSP3 | `Descriptors.FractionCSP3` | 0-1 | >=0.25 (Lovering 2009 escape-from-flatland) |
| QED | `QED.qed` | 0-1 | >=0.5 generally drug-like |
| SAscore | `sascorer.calculateScore` (external) | 1-10 | <=4 acceptable; >6 hard to synth |

**Goal:** Compute a standard physicochemical descriptor panel for drug-likeness filtering and QSAR features.

**Approach:** Combine RDKit `Descriptors`, `Lipinski`, and `QED` calls into a single dict so the caller gets MW, LogP, HBD/HBA, TPSA, rotatable bonds, aromatic rings, fraction sp3, and QED in one pass.

```python
from rdkit.Chem import Descriptors, Lipinski, QED

def physchem(mol):
    return {
        'MolWt': Descriptors.MolWt(mol),
        'MolLogP': Descriptors.MolLogP(mol),
        'HBD': Lipinski.NumHDonors(mol),
        'HBA': Lipinski.NumHAcceptors(mol),
        'TPSA': Descriptors.TPSA(mol),
        'RotBonds': Lipinski.NumRotatableBonds(mol),
        'AromRings': Lipinski.NumAromaticRings(mol),
        'FractionCSP3': Descriptors.FractionCSP3(mol),
        'QED': QED.qed(mol),
    }
```

## Drug-Likeness Rule Sets

| Rule | Constraints | Source |
|------|-------------|--------|
| Lipinski Ro5 | MW<=500, LogP<=5, HBD<=5, HBA<=10 | Lipinski 1997 |
| Veber | RotBonds<=10, TPSA<=140 | Veber 2002 (oral) |
| Ghose | 160<=MW<=480, -0.4<=LogP<=5.6, 40<=MR<=130, 20<=atoms<=70 | Ghose 1999 |
| Egan | LogP<=5.88, TPSA<=131.6 | Egan 2000 |
| Muegge | 200<=MW<=600, -2<=LogP<=5, TPSA<=150, rings<=7 | Muegge 2001 |
| Lead-like | MW<=350, LogP<=3 | Teague 1999 |
| Fragment Ro3 | MW<=300, LogP<=3, HBD<=3, HBA<=3, RotBonds<=3 | Congreve 2003 |
| BBB+ Pfizer CNS | TPSA<=90, MW<=500, HBD<=3 | Wager 2010 |

**Use case:** Apply Ro5/Veber as a screening filter, not a hard cutoff. ~30% of marketed oral drugs violate at least one Ro5 rule (analyzed by Doak 2014). For oncology indications, Ro5 deviation is common and acceptable.

## QED (Weighted Drug-Likeness)

QED (Bickerton 2012) is a single-number drug-likeness measure (0-1) combining 8 properties (MW, LogP, HBD, HBA, PSA, RotBonds, AromaticRings, structural alerts) via desirability functions.

**Caveat:** QED is trained on FDA-approved drugs; it under-rates fragment-like and natural-product-like molecules. Do not use as the sole drug-likeness filter for fragment screens or natural-product libraries.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Fingerprint changes between runs | Random seed not set for canonicalization | RDKit Morgan is deterministic; check if input differs (stereo, charges) |
| MACCS bit count != 166 | RDKit MACCS returns 167 bits (bit 0 unused) | Slice `[1:]` if comparing to literature 166-bit |
| Crippen LogP differs from XLogP | Different model | Use `Descriptors.MolLogP` for Crippen; `XLogP3` requires external lib |
| 3D descriptor differs between calls | Different conformer | Set `confId=0` explicitly; or average over ensemble |
| QED returns nan | Charged species or non-standard atom | Standardize (uncharge) before QED |
| Tanimoto on count vector wrong | `TanimotoSimilarity` expects bit vector | Use Hamming or weighted Tanimoto for counts |
| MolWt off by ~1 from PubChem | Implicit H counted differently | Use `Descriptors.ExactMolWt` for monoisotopic; PubChem reports average |

## References

- Rogers & Hahn, *J. Chem. Inf. Model.* 50:742 (2010) -- ECFP / Morgan fingerprints.
- Probst & Reymond, *J. Cheminformatics* 10:66 (2018) -- MHFP6 fingerprint.
- Capecchi et al., *J. Cheminformatics* 12:43 (2020) -- MAP4 fingerprint.
- Bickerton et al., *Nat. Chem.* 4:90 (2012) -- QED weighted drug-likeness.
- Lipinski et al., *Adv. Drug Deliv. Rev.* 23:3 (1997) -- Rule of 5.
- Veber et al., *J. Med. Chem.* 45:2615 (2002) -- Oral bioavailability rules.
- Lovering et al., *J. Med. Chem.* 52:6752 (2009) -- Fraction sp3 / escape from flatland.

## Related Skills

- chemoinformatics/molecular-io - Parse molecules before featurization
- chemoinformatics/molecular-standardization - Canonicalize before fingerprinting
- chemoinformatics/conformer-generation - Generate 3D for conformer-dependent descriptors
- chemoinformatics/similarity-searching - Use fingerprints for similarity ranking
- chemoinformatics/qsar-modeling - ML using these descriptors as features
- chemoinformatics/admet-prediction - Filter by drug-likeness criteria
- machine-learning/biomarker-discovery - ML on molecular features
