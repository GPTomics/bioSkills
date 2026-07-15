---
name: bio-molecular-descriptors
description: Calculates molecular fingerprints (ECFP/Morgan, FCFP, MACCS, RDKit, AtomPair, TopologicalTorsion, Avalon, MAP4, MHFP6) and physicochemical descriptors (Lipinski, QED, TPSA, Crippen LogP, 3D shape) with explicit choice tables, bit vs count semantics, and partial-charge model selection. Use when featurizing molecules for similarity, QSAR, virtual screening, or ML, or selecting the correct fingerprint for a chemotype-aware task.
tool_type: python
primary_tool: RDKit
---

## Version Compatibility

Reference examples tested with: RDKit 2024.09+, numpy 1.26+, pandas 2.2+, map4 1.1+ (MAP4), mhfp 1.9+. Use `mapchiral` separately when the stereochemistry-aware MAP4C fingerprint is intended.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Molecular Descriptors

Featurize molecules for similarity search, QSAR, virtual screening, or ML. Fingerprint performance is **dataset- and objective-dependent**: ECFP4 is a strong drug-like baseline, atom-pair and topological-torsion fingerprints expose longer-range topology, MAP4/MHFP6 target broader chemical-space searches, and 3D conformer-based descriptors are needed when shape and stereochemistry matter.

For canonicalization before featurization, see `chemoinformatics/molecular-standardization`. For 3D-only descriptors, see `chemoinformatics/conformer-generation`.

## Fingerprint Taxonomy

| Fingerprint | Type | Radius/Path | Bits | Use case | Fails when |
|-------------|------|-------------|------|----------|------------|
| Morgan (ECFP) | Circular | r=2 (ECFP4), r=3 (ECFP6) | 2048 typical | Drug-like similarity, ML default | Loses long-range topology; bit collisions at low nBits |
| FCFP | Functional Morgan | r=2 default | 2048 | Pharmacophore-aware similarity | Same caveats as ECFP; less specific |
| MACCS | Substructure key | 166 fixed bits | 167 | Quick fingerprint, drug-likeness | Too sparse for large diverse libraries |
| RDKit FP | Path/subgraph-based | paths and branched subgraphs up to 7 bonds by default | 2048 | RDKit-native ECFP alternative | Drug-like only; not optimal for scaffold hopping |
| AtomPair | Pair + topological distance | All atom pairs | 2048 | Long-range topological similarity | Slower than ECFP; harder to interpret |
| TopologicalTorsion | 4-atom torsion | All TT | 2048 | Path-pattern similarity | Like AP, slower than ECFP |
| Avalon | Substructure + atom pairs | Mixed | 512/1024 | Fast similarity | Less standard; older |
| MAP4 (MinHashed atom-pair) | MinHash atom-pair | r=1,2 | 1024/2048 | Biological + metabolite diversity | `map4` library required; slower hash |
| MHFP6 (MinHash) | MinHash ECFP-like | r=3 (diam 6) | 2048 | Large-library nearest-neighbor with a compatible MinHash/LSH index | Different distance semantics from folded-bit Tanimoto |
| Pharm2D | 2D pharmacophore | feature pairs/triplets | sparse | Pharmacophore search | Sparse, slower |

**Decision:** For drug-like similarity ranking, start with **ECFP4 2048 bit** because it is fast and well characterized. MHFP6 outperformed ECFP4 for analog recovery in the benchmark reported by Probst and Reymond (2018), making it a candidate for large, diverse libraries. For scaffold hopping, benchmark ECFP4, AtomPair, TopologicalTorsion, and pharmacophore fingerprints on target-relevant actives and decoys; published comparisons do not support a universal AtomPair advantage (Gardiner et al. 2011; Riniker & Landrum 2013).

## Bit vs Count Vectors

| Form | Use | Library impact |
|------|-----|----------------|
| Bit (0/1) | Tanimoto similarity, BulkTanimotoSimilarity, RDKit fingerprint folding | Standard for similarity |
| Count (integer) | Some ML methods, RF on counts, neural fingerprints | Loses bit-level fast operations; richer signal |
| Sparse (dict) | Direct chemical interpretation (which fragments at which atoms) | Use for SHAP / atomic attribution |

```python
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

mol = Chem.MolFromSmiles('CCO')

morgan = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
ecfp4_bit = morgan.GetFingerprint(mol)
ecfp4_count = morgan.GetCountFingerprint(mol)
ecfp4_sparse = morgan.GetSparseCountFingerprint(mol)
```

## Morgan / ECFP Radius Math

ECFP-X notation: X is the **diameter** in bonds. RDKit's `radius` parameter is half of X.

| Notation | RDKit radius | Diameter | Captures |
|----------|--------------|----------|----------|
| ECFP0 | 0 | 0 | Atom identity only |
| ECFP2 | 1 | 2 | Atom + immediate neighbors |
| ECFP4 | 2 | 4 | Atom + 2-bond environment |
| ECFP6 | 3 | 6 | Atom + 3-bond environment |

**Trade-off:** Larger radius captures more specific local environments but increases collisions at fixed `nBits`. **ECFP4 2048** is a common baseline (Rogers & Hahn 2010; Wu et al. 2018). O'Boyle and Sayle (2016) showed that increasing folded-fingerprint length can improve virtual-screening performance, but they did not establish a universal 4096-bit setting or a 1-5% collision rate. Measure collision occupancy and model performance for the dataset; increase `nBits` or use an unhashed sparse representation when needed.

## FCFP vs ECFP

FCFP (Functional-Class) uses RDKit's Morgan feature invariants (donor, acceptor, aromatic, halogen, basic, and acidic) instead of atom identity. Hydrophobe is a family in `BaseFeatures.fdef`, but it is not one of the default Morgan feature-invariant classes. FCFP trades atom-specificity for functional-equivalence.

```python
ecfp_generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
feature_invariants = rdFingerprintGenerator.GetMorganFeatureAtomInvGen()
fcfp_generator = rdFingerprintGenerator.GetMorganGenerator(
    radius=2, fpSize=2048, atomInvariantsGenerator=feature_invariants)
ecfp4 = ecfp_generator.GetFingerprint(mol)
fcfp4 = fcfp_generator.GetFingerprint(mol)
```

**When to use FCFP4:** Scaffold-hopping campaigns, pharmacophore-driven similarity, cross-target activity prediction.

**When to use ECFP4:** Within-series QSAR, lead optimization, when chemotype identity matters.

## 3D Descriptors and Conformer Dependence

Conformer-dependent descriptors (asphericity, eccentricity, principal moments of inertia, RDF) require a generated 3D structure. A single conformer may be unrepresentative when the molecule is flexible; measure descriptor variation across a conformer ensemble when the downstream conclusion depends on 3D shape.

**Goal:** Compute 3D shape descriptors over a conformer ensemble rather than from a single (possibly unrepresentative) conformer.

**Approach:** Add explicit hydrogens, embed N conformers with ETKDGv3, MMFF-optimize them all, then evaluate the descriptor across each conformer for downstream averaging.

```python
from rdkit.Chem import AllChem, Descriptors3D

mol = Chem.MolFromSmiles('CCCCO')
mol = Chem.AddHs(mol)

params = AllChem.ETKDGv3()
params.randomSeed = 42
conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=20, params=params)
if not conf_ids:
    raise RuntimeError('ETKDGv3 failed to generate any conformers')
if not AllChem.MMFFHasAllMoleculeParams(mol):
    raise ValueError('MMFF94 parameters are unavailable for this molecule')
optimization_results = AllChem.MMFFOptimizeMoleculeConfs(mol)
if any(status != 0 for status, _ in optimization_results):
    raise RuntimeError('MMFF94 optimization did not converge for every conformer')

asphericities = [Descriptors3D.Asphericity(mol, confId=c) for c in conf_ids]
```

**Decision:** For QSAR / ML, choose and document the conformer count using a convergence check on representative molecules. Report the aggregation rule, such as a simple mean or a Boltzmann-weighted average, and the energy model used for any weights.

## Partial Charge Methods

| Method | Software | Cost | Accuracy | Use for |
|--------|----------|------|----------|---------|
| Gasteiger-Marsili | RDKit, Open Babel | Fast | Empirical, rough | Charge-aware preparation or models that explicitly require Gasteiger charges; Vina/Vinardo scoring itself does not require assigned atom charges |
| MMFF94 | RDKit | 0.1s/mol | Force-field consistent | MMFF energy, conformer ranking |
| AM1-BCC | antechamber (AmberTools) | ~10s/mol | Semi-empirical | MD setup, FEP, GAFF |
| RESP | psi4, Gaussian | minutes/mol | Restrained fit to a quantum-mechanical ESP; protocol-specific | Force-field workflows parameterized for that RESP protocol |
| OpenFF Recharge | openff-recharge | Workflow-dependent | Framework for generating/retrieving QC ESP data and fitting library charges, BCCs, RESP charges, or virtual sites | Developing or evaluating charge models; it is not one charge-assignment method |

```python
from rdkit.Chem import AllChem

AllChem.ComputeGasteigerCharges(mol)
for atom in mol.GetAtoms():
    print(atom.GetIdx(), atom.GetPropsAsDict().get('_GasteigerCharge', None))
```

**Critical:** Charge method must match downstream. Gasteiger charges in an AMBER MD run violate the assumptions of the protein force field.

## MAP4 and MHFP6 for Diverse Libraries

For libraries spanning drug-like molecules, natural products, peptides, and metabolites, compare ECFP4 with MAP4 or MHFP6 on task-relevant retrieval benchmarks. MAP4 and MHFP6 use MinHash with atom-pair or circular-substructure shingles, but no universal pairwise-similarity range establishes that ECFP4 is saturated for every mixed library.

```python
from mhfp.encoder import MHFPEncoder

encoder = MHFPEncoder(2048)
mhfp6 = encoder.encode_mol(mol, radius=3)
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
| FractionCSP3 | `Descriptors.FractionCSP3` | 0-1 | Descriptive; higher sp3 character was associated with clinical progression by Lovering et al. (2009), without a universal cutoff |
| QED | `QED.qed` | 0-1 | Higher is more similar to the reference property distributions; a project may use >=0.5 as a triage heuristic |
| SAscore | `sascorer.calculateScore` (external) | 1-10 | Lower is easier by the model; project cutoffs such as <=4 or >6 require dataset calibration |

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
| Muegge | 200<=MW<=600, -2<=LogP<=5, TPSA<=150, rings<=7, C>4, heteroatoms>1, RotBonds<=15, HBD<=5, HBA<=10 | Muegge 2001 |
| Lead-like | MW<=350, LogP<=3 | Teague 1999 |
| Fragment Ro3 | MW<=300, LogP<=3, HBD<=3, HBA<=3, RotBonds<=3, TPSA<=60 A^2 | Congreve 2003 |
| Pfizer CNS MPO | Six desirability functions: ClogP, ClogD, MW, TPSA, HBD, and pKa | Wager 2010 |

**Use case:** Treat Ro5 and Veber criteria as risk indicators rather than universal hard cutoffs. Doak et al. (2014) analyze orally bioavailable drugs and candidates beyond the Rule of 5, but do not support the claim that approximately 30% of marketed oral drugs violate at least one rule. For CNS prioritization, implement the six-property Wager MPO desirability score rather than replacing it with three hard thresholds.

## QED (Weighted Drug-Likeness)

QED (Bickerton 2012) is a single-number drug-likeness measure (0-1) combining 8 properties (MW, LogP, HBD, HBA, PSA, RotBonds, AromaticRings, structural alerts) via desirability functions.

**Caveat:** QED summarizes desirability functions derived from property distributions of marketed oral drugs; it is not a supervised predictor trained specifically on FDA-approved drugs. It can under-rank fragment-like or natural-product-like molecules, so do not use it as the sole filter for those libraries.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Fingerprint changes between runs | Random seed not set for canonicalization | RDKit Morgan is deterministic; check if input differs (stereo, charges) |
| MACCS bit count != 166 | RDKit MACCS returns 167 bits (bit 0 unused) | Slice `[1:]` if comparing to literature 166-bit |
| Crippen LogP differs from XLogP | Different model | Use `Descriptors.MolLogP` for Crippen; `XLogP3` requires external lib |
| 3D descriptor differs between calls | Different conformer | Set `confId=0` explicitly; or average over ensemble |
| QED returns nan | Charged species or non-standard atom | Standardize (uncharge) before QED |
| Count-vector similarity differs from bit-vector similarity | Count multiplicities change the generalized Tanimoto calculation | RDKit supports Tanimoto on sparse count vectors; record the vector type and do not compare its threshold directly with a folded-bit threshold |
| MolWt off by ~1 from PubChem | Implicit H counted differently | Use `Descriptors.ExactMolWt` for monoisotopic; PubChem reports average |

## References

- Rogers & Hahn, *J. Chem. Inf. Model.* 50:742-754 (2010) -- ECFP fingerprints. https://doi.org/10.1021/ci100050t
- Probst & Reymond, *J. Cheminformatics* 10:66 (2018) -- MHFP6 fingerprint. https://doi.org/10.1186/s13321-018-0321-8
- Capecchi et al., *J. Cheminformatics* 12:43 (2020) -- MAP4 fingerprint. https://doi.org/10.1186/s13321-020-00445-4
- MAP4, official package -- installation and current Python interface. https://pypi.org/project/map4/
- OpenFF Recharge, official documentation -- supported charge-model generation and fitting workflows. https://docs.openforcefield.org/projects/recharge/en/stable/
- AutoDock Vina, official documentation -- Vina/Vinardo charge semantics. https://autodock-vina.readthedocs.io/en/stable/
- Gardiner et al., *Future Med. Chem.* 3:405-414 (2011) -- scaffold-hopping fingerprint comparison. https://doi.org/10.4155/fmc.11.4
- Riniker & Landrum, *J. Cheminformatics* 5:26 (2013) -- fingerprint benchmarking. https://doi.org/10.1186/1758-2946-5-26
- O'Boyle & Sayle, *J. Cheminformatics* 8:36 (2016) -- fingerprint folding and virtual-screening performance. https://doi.org/10.1186/s13321-016-0148-0
- Wu et al., *Chem. Sci.* 9:513-530 (2018) -- MoleculeNet benchmarks. https://doi.org/10.1039/C7SC02664A
- Bickerton et al., *Nat. Chem.* 4:90-98 (2012) -- QED weighted drug-likeness. https://doi.org/10.1038/nchem.1243
- Lipinski et al., *Adv. Drug Deliv. Rev.* 23:3-25 (1997) -- Rule of 5. https://doi.org/10.1016/S0169-409X(96)00423-1
- Veber et al., *J. Med. Chem.* 45:2615-2623 (2002) -- oral bioavailability criteria. https://doi.org/10.1021/jm020017n
- Ghose et al., *J. Comb. Chem.* 1:55-68 (1999) -- physicochemical property ranges. https://doi.org/10.1021/cc9800071
- Egan et al., *J. Med. Chem.* 43:3867-3877 (2000) -- absorption model. https://doi.org/10.1021/jm000292e
- Muegge et al., *J. Med. Chem.* 44:1841-1846 (2001) -- drug-like property filters. https://doi.org/10.1021/jm015507e
- Teague et al., *Angew. Chem. Int. Ed.* 38:3743-3748 (1999) -- lead-like libraries. https://doi.org/10.1002/%28SICI%291521-3773%2819991216%2938%3A24%3C3743%3A%3AAID-ANIE3743%3E3.0.CO%3B2-U
- Congreve et al., *Drug Discov. Today* 8:876-877 (2003) -- Rule of 3. https://doi.org/10.1016/S1359-6446(03)02831-9
- Lovering et al., *J. Med. Chem.* 52:6752-6756 (2009) -- fraction sp3 and clinical progression. https://doi.org/10.1021/jm901241e
- Ritchie & Macdonald, *Drug Discov. Today* 14:1011-1020 (2009) -- aromatic ring count and developability. https://doi.org/10.1016/j.drudis.2009.07.014
- Ertl & Schuffenhauer, *J. Cheminformatics* 1:8 (2009) -- synthetic accessibility score. https://doi.org/10.1186/1758-2946-1-8
- Wager et al., *ACS Chem. Neurosci.* 1:435-449 (2010) -- CNS multiparameter optimization. https://doi.org/10.1021/cn100008c
- Doak et al., *Chem. Biol.* 21:1115-1142 (2014) -- orally bioavailable drugs beyond the Rule of 5. https://doi.org/10.1016/j.chembiol.2014.08.013

## Related Skills

- chemoinformatics/molecular-io - Parse molecules before featurization
- chemoinformatics/molecular-standardization - Canonicalize before fingerprinting
- chemoinformatics/conformer-generation - Generate 3D for conformer-dependent descriptors
- chemoinformatics/similarity-searching - Use fingerprints for similarity ranking
- chemoinformatics/qsar-modeling - ML using these descriptors as features
- chemoinformatics/admet-prediction - Filter by drug-likeness criteria
- machine-learning/biomarker-discovery - ML on molecular features
