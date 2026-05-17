---
name: bio-pharmacophore-modeling
description: Builds and applies 3D pharmacophore models using RDKit Pharm3D, apo2ph4 (receptor-based), Pharmer / Pharmit (search), and PharmacoForge (diffusion-based generation), covering ligand-based pharmacophore (from active set alignment) and receptor-based pharmacophore (from binding pocket geometry). Explicit handling of feature types, geometric tolerances, partial matching, and pharmacophore-based virtual screening. Use when identifying scaffold-hopping candidates, building shape-and-feature search queries, or transferring SAR across chemotypes.
tool_type: python
primary_tool: RDKit
---

## Version Compatibility

Reference examples tested with: RDKit 2024.09+, pharmer / pharmit (web service), PharmIT 1.1+, plip 2.4+ (interaction analysis).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show rdkit` then `help(rdkit.Chem.Pharm3D)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Pharmacophore Modeling

Build 3D pharmacophore queries that capture the essential interaction features of a ligand-target binding event. A pharmacophore is the *spatial arrangement of pharmacophore features* (donor, acceptor, hydrophobe, aromatic, charged) sufficient for activity, abstracted from any specific chemotype. Used for scaffold-hopping (find compounds with different scaffold but matching pharmacophore), virtual screening (faster than docking), and cross-target SAR transfer. Modern best practice: derive pharmacophore from co-crystal structure if available (receptor-based with apo2ph4) or align actives if no crystal (ligand-based). Diffusion-based generation (PharmacoForge, Karki 2025) lets pharmacophore drive de novo design.

For 2D scaffold-based searches, see `chemoinformatics/scaffold-analysis`. For 3D shape similarity, see `chemoinformatics/shape-similarity`. For protein-ligand interaction analysis, see `chemoinformatics/virtual-screening`.

## Pharmacophore Feature Types

| Feature | RDKit code | Definition | Geometric tolerance |
|---------|------------|------------|----------------------|
| H-bond donor | D | -OH, -NH | 1.0-1.5 Å |
| H-bond acceptor | A | sp2 O / N (lone pair) | 1.0-1.5 Å |
| Hydrophobe | H | sp3 C / aromatic ring centroid | 1.5-2.0 Å |
| Aromatic ring | R | Aromatic ring centroid + normal | 1.0-1.5 Å |
| Positive ionizable | P | -NH3+, -NR3+ | 1.0-1.5 Å |
| Negative ionizable | N | -COO-, -SO3- | 1.0-1.5 Å |
| Halogen | X | Cl, Br, I (halogen bond donor) | 1.0-1.5 Å |
| Metal coordination | M | sp/sp2 N/O near metal | 0.5-1.0 Å |

Tolerances are pharmacophore-feature distance windows in the search. Tighter tolerances = fewer hits but more specific.

## Method Taxonomy

| Method | Origin | Use case | Fails when |
|--------|--------|----------|------------|
| Ligand-based (LBP) | Catalyst, MOE, RDKit Pharm3D | Multiple actives, no crystal | <3 actives; flexible actives |
| Receptor-based (RBP) | apo2ph4, LigandScout | Co-crystal available | Apo structure (use AlphaFold3 or Boltz) |
| Common pharmacophore | Pharm3D `EmbedPharmacophore` | Consensus from active set | Diverse actives confound alignment |
| Diffusion-based (PharmacoForge) | Karki 2025 | De novo generation with pharmacophore prior | Pretrained model required |
| Active learning pharmacophore | Catalyst variant | Iterative refinement | Custom; not standard |

## Decision Tree by Scenario

| Scenario | Method | Tools |
|----------|--------|-------|
| Co-crystal structure available | Receptor-based | apo2ph4 + Pharmer/Pharmit |
| Multiple active compounds, no crystal | Ligand-based common pharmacophore | RDKit Pharm3D `EmbedPharmacophore` |
| Single active compound | Single-conformer pharmacophore | RDKit Pharm3D from bioactive conformer |
| Scaffold hopping prospective | Receptor-based + shape filter | apo2ph4 + ROCS |
| Cross-target SAR transfer | Common pharmacophore across targets | Manual + LigandScout |
| De novo design with pharmacophore | PharmacoForge | Diffusion-based generation |
| Library pre-filtering | Pharmacophore screen | Pharmit search |

## Ligand-Based Pharmacophore (RDKit Pharm3D)

**Goal:** Extract a common pharmacophore from a set of bioactive compounds.

**Approach:** Align actives by maximum-common-substructure or shape; extract conserved features; derive a pharmacophore signature.

```python
from rdkit import Chem
from rdkit.Chem import AllChem, ChemicalFeatures
from rdkit.Chem.Pharm3D import Pharmacophore
from rdkit.RDPaths import RDDataDir
import os

fdef_file = os.path.join(RDDataDir, 'BaseFeatures.fdef')
factory = ChemicalFeatures.BuildFeatureFactory(fdef_file)

active_smiles = ['CC(C)c1ccc(C(=O)NCc2ccccn2)cc1', 'CCC(C)c1ccc(C(=O)NCc2ccccn2)cc1']
active_mols = [Chem.AddHs(Chem.MolFromSmiles(s)) for s in active_smiles]
for m in active_mols:
    AllChem.EmbedMolecule(m, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(m)

# Extract pharmacophore features (family/type/atom-ids/3D-pos) per molecule
feature_lists = [factory.GetFeaturesForMol(m) for m in active_mols]

# Build a target pharmacophore from one active's features; the bounds matrix
# encodes per-feature-pair distance ranges across all conformer/active variability.
# Below uses 2 features (Aromatic + Donor) with 2.5-3.5 A and 5.0-6.0 A bands.
feature_types = ['Aromatic', 'Donor']
# bounds is symmetric 2x2 distance-range matrix; off-diagonal = (lo, hi)
import numpy as np
bounds_matrix = np.array([[0.0, 5.0], [3.5, 0.0]])  # upper / lower bounds
pharmacophore = Pharmacophore.Pharmacophore(feature_types)
pharmacophore.setLowerBound(0, 1, 3.5)
pharmacophore.setUpperBound(0, 1, 5.0)

# To search a target molecule for matches, use Pharm3D.EmbedLib.EmbedPharmacophore
# (see chemoinformatics/conformer-generation for 3D embedding fundamentals)
```

`BaseFeatures.fdef` (RDKit-shipped) defines feature SMARTS. For drug-like pharmacophores, this is the standard starting point.

## Receptor-Based Pharmacophore (apo2ph4)

**Goal:** Derive a pharmacophore from a protein binding-pocket structure without requiring a bound ligand.

**Approach:** Identify hot-spots (donor / acceptor / hydrophobe regions) from protein geometry; assemble into a pharmacophore.

```bash
apo2ph4 -pdb receptor.pdb -site_residues 'A:100,A:101,A:104,A:108' \
        -output pharmacophore.ph4
```

apo2ph4 outputs a `.ph4` format compatible with Phase, MOE, and Pharmer.

When a co-crystal ligand is available, **derive pharmacophore directly from the ligand binding pose**: each ligand feature in contact with a complementary protein residue is part of the pharmacophore.

```python
from plip.basic import config
from plip.structure.preparation import PDBComplex

mol_complex = PDBComplex()
mol_complex.load_pdb('complex.pdb')
mol_complex.analyze()

for site in mol_complex.interaction_sets.values():
    for interaction in site.all_itypes:
        # H-bond donor / acceptor / pi-stacking / hydrophobic / salt-bridge
        feature_type = interaction.type
        feature_atom_coords = interaction.ligatom.coords
```

PLIP outputs interaction types per ligand atom; each interaction maps to a pharmacophore feature.

## Pharmacophore Search (Pharmit / Pharmer)

For library screening, Pharmer/Pharmit are the standard tools.

```bash
# Pharmer command line (offline alternative)
pharmer search -q pharmacophore.ph4 -dbdir zinc_db -out hits.sdf
```

Or use Pharmit web service (https://pharmit.csb.pitt.edu) for browser-based or REST search.

```python
import requests

url = 'https://pharmit.csb.pitt.edu/...'  # specific endpoint
response = requests.post(url, json={'pharmacophore': pharmacophore_dict})
hits = response.json()
```

Pharmacophore screen is orders of magnitude faster than docking; typical 100M-compound search in minutes.

## Pharmacophore Quality Validation

Evaluate a pharmacophore by:

1. **Retrospective enrichment**: AUC-like metric on actives vs decoys (DUD-E, COCONUT)
2. **Geometric tightness**: feature distance variance across actives
3. **Selectivity**: false positives in inactive set should be low
4. **Specific consistency**: pharmacophore matches each active's bioactive conformer

```python
def pharmacophore_enrichment(query_pharmacophore, actives, inactives):
    n_active_match = sum(matches_pharmacophore(a, query_pharmacophore) for a in actives)
    n_inactive_match = sum(matches_pharmacophore(d, query_pharmacophore) for d in inactives)
    enrichment = (n_active_match / len(actives)) / (n_inactive_match / len(inactives))
    return enrichment
```

A good pharmacophore: enrichment >= 5x relative to random.

## Diffusion-Based Pharmacophore Design (PharmacoForge)

PharmacoForge (Karki 2025) generates molecules conditioned on a pharmacophore query using a diffusion model:

```python
# Pseudo-code; depends on PharmacoForge installation
# from pharmacoforge import PharmacophoreDiffusionGenerator
# gen = PharmacophoreDiffusionGenerator()
# generated = gen.generate(pharmacophore_query, n_molecules=100)
```

Trade-off: PharmacoForge produces de novo molecules that satisfy pharmacophore; lower drug-likeness than REINVENT but better novelty.

## Pharmacophore vs Shape vs 2D Fingerprint

| Method | Captures | Best for |
|--------|----------|----------|
| ECFP4 Tanimoto | Local atom environments | Lead optimization (same series) |
| FCFP4 Tanimoto | Pharmacophore-equivalent atoms | Loose similarity in series |
| Shape similarity (ROCS) | 3D shape volume | Scaffold hopping by shape |
| Pharmacophore | Discrete features in space | Scaffold hopping with feature specificity |
| Combined (Tanimoto + shape) | Multi-objective | Production VS |

Pharmacophore is more *interpretable* than shape: a hit explains why it matched (donor at position X, hydrophobe at position Y).

## Per-Tool Failure Modes

### Ligand-based -- diverse actives confound

**Trigger:** Active set spans multiple scaffolds with different bound conformations.

**Mechanism:** No common pharmacophore exists; algorithm forces non-consensus features.

**Symptom:** Pharmacophore matches no actives in retrospective.

**Fix:** Cluster actives by scaffold first; derive per-cluster pharmacophore.

### Receptor-based -- apo structure

**Trigger:** Protein in apo form (no bound ligand).

**Mechanism:** Side-chain rotamers differ between apo and holo; "binding site" geometry is wrong.

**Symptom:** Pharmacophore inferred from apo doesn't match holo experimental data.

**Fix:** Use AlphaFold3 / Boltz-1 to predict holo conformation; derive pharmacophore from predicted holo.

### Pharmacophore -- single conformer bias

**Trigger:** Active aligned to its first generated conformer, not bioactive conformer.

**Mechanism:** Crystal structure not available; generated conformer may not be the bound one.

**Symptom:** Pharmacophore inconsistent across runs (different starting conformer chosen).

**Fix:** Use conformer ensemble; align all to common scaffold; choose conformer most consistent with other actives.

### Tolerance too tight

**Trigger:** Default geometric tolerance < 0.5 Å.

**Mechanism:** Real bioactive conformers have flexibility; rigid pharmacophore filters most molecules out.

**Symptom:** Search returns zero hits.

**Fix:** Use tolerance 1.0-1.5 Å for drug-like; up to 2 Å for flexible peptide-like.

### Pharmacophore search misses bioisostere

**Trigger:** Bioisostere replacement (e.g., -COOH replaced by tetrazole).

**Mechanism:** Tetrazole functions as acid bioisostere but RDKit features may not classify identically.

**Symptom:** Known bioisosteric active not found.

**Fix:** Use ChemAxon-style bioisosteric feature equivalence; or pharmacophore feature class expansion (acid generic vs -COOH specific).

### PLIP -- water-mediated interaction missed

**Trigger:** Bridging water between ligand donor and protein acceptor.

**Mechanism:** PLIP doesn't include explicit water.

**Symptom:** Pharmacophore missing critical H-bond feature.

**Fix:** Add water-mediated interaction manually based on crystal water positions; or use PoseView for full interaction view.

## Reconciliation: Ligand-Based vs Receptor-Based

| Aspect | Ligand-based | Receptor-based |
|--------|--------------|----------------|
| Data needed | Multiple actives | Co-crystal or predicted holo |
| Bias | Toward known chemotype | None |
| Hit set | Similar to known actives | More diverse |
| Discovery potential | Limited | Higher |
| Pharmacophore confidence | Higher (validated against multiple actives) | Lower (single ligand) |

For prospective scaffold-hopping, receptor-based is preferred. For ranking analogs, ligand-based is sufficient.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `Pharm3D.EmbedPharmacophore` fails | Bounds matrix infeasible | Loosen tolerances; increase max_attempts |
| Pharmacophore matches everything | Too few features | Add features; tighten tolerances |
| Pharmacophore matches nothing | Too many features or tight bounds | Reduce feature count; loosen tolerances |
| BaseFeatures.fdef not found | RDKit installation issue | Check `from rdkit.RDPaths import RDDataDir` |
| Pharmacophore-conformer mismatch | Wrong conformer used | Use bioactive conformer from crystal |
| Pharmit search timeout | Library too large | Pre-filter by 2D fingerprint Tanimoto |
| apo2ph4 produces empty | No druggable hot-spots | Lower thresholds; use shape-only filter |

## References

- Wolber & Langer, *J. Chem. Inf. Model.* 45:160 (2005) -- LigandScout pharmacophore.
- Sander et al., *J. Chem. Inf. Model.* (2024) -- apo2ph4 for receptor-based.
- Karki et al., *Front. Bioinform.* (2025) -- PharmacoForge diffusion pharmacophore.
- Stiefl et al., *J. Chem. Inf. Model.* (2021) -- RDKit pharmacophore framework.
- Adasme et al., *Nucleic Acids Res.* 49:W530 (2021) -- PLIP interaction profiler.
- Koes et al., *Nucleic Acids Res.* 44:W436 (2016) -- Pharmit interactive search.

## Related Skills

- chemoinformatics/molecular-io - Parse molecules
- chemoinformatics/conformer-generation - Generate 3D for pharmacophore
- chemoinformatics/shape-similarity - 3D shape adjacent to pharmacophore
- chemoinformatics/virtual-screening - Pharmacophore as docking pre-filter
- chemoinformatics/scaffold-analysis - 2D scaffold-hopping context
- chemoinformatics/generative-design - PharmacoForge for de novo
- structural-biology/structure-io - PDB handling
