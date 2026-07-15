---
name: bio-pharmacophore-modeling
description: Builds and applies 3D pharmacophore models using RDKit Pharm3D, the apo2ph4 receptor-based workflow (Heider et al. 2023), Pharmer / Pharmit for search, and PharmacoForge for protein-pocket-conditioned pharmacophore generation (Flynn et al. 2025), covering ligand-based pharmacophores from active-set alignment and receptor-based pharmacophores from binding-pocket geometry. Explicitly handles feature types, geometric tolerances, partial matching, and pharmacophore-based virtual screening. Use when identifying scaffold-hopping candidates, building shape-and-feature search queries, or transferring SAR across chemotypes.
tool_type: python
primary_tool: RDKit
---

## Version Compatibility

Reference examples tested with: RDKit 2024.09+, Pharmit web service, and PLIP 2.4+ (interaction analysis). Verify the deployed Pharmit/Pharmer interface and query format before automation.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show rdkit` then `help(rdkit.Chem.Pharm3D)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Pharmacophore Modeling

Build 3D pharmacophore queries that capture the essential interaction features of a ligand-target binding event. A pharmacophore is the *spatial arrangement of pharmacophore features* (donor, acceptor, hydrophobe, aromatic, charged) sufficient for activity, abstracted from any specific chemotype. Use pharmacophores for scaffold hopping, virtual-screening prefilters, and cross-target SAR transfer. Derive interaction features directly from a co-crystal when available, use apo2ph4 to derive models from an apo pocket (Heider et al. 2023), or align known actives for a ligand-based model. PharmacoForge generates candidate 3D pharmacophores conditioned on a protein pocket; those pharmacophores can then retrieve matching molecules from a library (Flynn et al. 2025).

For 2D scaffold-based searches, see `chemoinformatics/scaffold-analysis`. For 3D shape similarity, see `chemoinformatics/shape-similarity`. For protein-ligand interaction analysis, see `chemoinformatics/virtual-screening`.

## Pharmacophore Feature Types

| Feature | Common shorthand | Definition | Geometric tolerance |
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

The ranges in this table are repository starting heuristics, not universal feature tolerances. Set final bounds from aligned-feature variability, coordinate uncertainty, and retrospective validation for the selected search engine.

The one-letter labels above are human-readable shorthand, not RDKit API codes. RDKit's shipped `BaseFeatures.fdef` uses family names such as `Donor`, `Acceptor`, `Hydrophobe`, `Aromatic`, `PosIonizable`, and `NegIonizable`. Its default feature definitions do not provide every halogen-bond or metal-coordination model; add and validate project-specific feature definitions when those interactions matter.

## Method Taxonomy

| Method | Origin | Use case | Fails when |
|--------|--------|----------|------------|
| Ligand-based (LBP) | Catalyst, MOE, RDKit Pharm3D | Multiple actives, no crystal | <3 actives; flexible actives |
| Receptor-based (RBP) | apo2ph4, LigandScout, PLIP | Co-crystal or a defined apo pocket | Uncertain pocket conformation |
| Common pharmacophore | Validated alignment/feature-consensus workflow; RDKit can represent and query the resulting model | Consensus from active set | Diverse actives or uncertain bioactive conformers confound alignment |
| Pocket-conditioned generation (PharmacoForge) | Flynn et al. 2025 | Generate candidate pharmacophores from a protein pocket | Does not directly generate molecules; pretrained model required |
| Active learning pharmacophore | Catalyst variant | Iterative refinement | Custom; not standard |

## Decision Tree by Scenario

| Scenario | Method | Tools |
|----------|--------|-------|
| Co-crystal structure available | Interaction-derived receptor model | PLIP or LigandScout + Pharmit |
| Apo structure with a defined pocket | Apo receptor model | apo2ph4; export LigandScout PML |
| Multiple active compounds, no crystal | Ligand-based common pharmacophore | Alignment plus consensus-feature derivation in validated custom or external tooling; RDKit Pharm3D can apply the resulting model |
| Single active compound | Single-conformer pharmacophore | RDKit Pharm3D from bioactive conformer |
| Scaffold hopping prospective | Receptor-based + shape filter | apo2ph4 or interaction-derived model + shape search |
| Cross-target SAR transfer | Common pharmacophore across targets | Manual + LigandScout |
| Generate pocket-conditioned pharmacophores | PharmacoForge | Diffusion model followed by library retrieval |
| Library pre-filtering | Pharmacophore screen | Pharmit search |

## Ligand-Based Pharmacophore (RDKit Pharm3D)

**Goal:** Derive a common pharmacophore from aligned bioactive conformers, then apply that established model to candidate molecules.

**Approach:** Consensus derivation is a separate modeling step: select or generate plausible bioactive conformers, align them using a documented method, identify conserved feature correspondences, and estimate distance bounds or tolerances. RDKit does not provide a single `EmbedPharmacophore` call that performs those steps. `EmbedPharmacophore` instead generates conformations of a molecule that satisfy an already defined pharmacophore.

```python
from rdkit import Chem, Geometry
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm3D import EmbedLib, Pharmacophore
from rdkit.RDPaths import RDDataDir
import os

fdef_file = os.path.join(RDDataDir, 'BaseFeatures.fdef')
factory = ChemicalFeatures.BuildFeatureFactory(fdef_file)

# This is an already defined model. Coordinates and bounds must come from a
# validated consensus-derivation workflow or another justified source. RDKit
# requires FreeChemicalFeature objects, not feature-family strings.
query_features = [
    ChemicalFeatures.FreeChemicalFeature(
        'Aromatic', Geometry.Point3D(0.0, 0.0, 0.0)),
    ChemicalFeatures.FreeChemicalFeature(
        'Donor', Geometry.Point3D(4.0, 0.0, 0.0)),
]
pharmacophore = Pharmacophore.Pharmacophore(query_features)
pharmacophore.setLowerBound(0, 1, 3.5)
pharmacophore.setUpperBound(0, 1, 5.0)

target = Chem.AddHs(Chem.MolFromSmiles('c1ccc(cc1)CCN'))
can_match, feature_matches = EmbedLib.MatchPharmacophoreToMol(
    target, factory, pharmacophore)
if can_match:
    atom_match = tuple(tuple(matches[0].GetAtomIds())
                       for matches in feature_matches)
    _, embeddings, n_failed = EmbedLib.EmbedPharmacophore(
        target, atom_match, pharmacophore, randomSeed=23, silent=True)
```

`BaseFeatures.fdef` (RDKit-shipped) defines feature SMARTS and is a useful starting feature taxonomy. The code above demonstrates applying an existing two-feature model; it does not infer a consensus model from active compounds.

## Receptor-Based Pharmacophore (apo2ph4 workflow)

**Goal:** Derive a pharmacophore from a protein binding-pocket structure without requiring a bound ligand.

**Approach:** Identify donor, acceptor, and hydrophobic hot spots from apo-pocket geometry, cluster them, and assemble candidate pharmacophores. Heider et al. describe apo2ph4 in *J. Chem. Inf. Model.* 63:101-110 (2023). Use the source release's documented scripts and environment rather than assuming a packaged `apo2ph4` command: the published workflow writes LigandScout PML output, not a generic `.ph4` file. Treat conversion to Pharmit, Pharmer, MOE, or Phase as a separate, explicitly validated step because pharmacophore formats are not interchangeable.

When a co-crystal ligand is available, **derive pharmacophore directly from the ligand binding pose**: each ligand feature in contact with a complementary protein residue is part of the pharmacophore.

```python
from plip.basic import config
from plip.structure.preparation import PDBComplex

mol_complex = PDBComplex()
mol_complex.load_pdb('complex.pdb')
mol_complex.analyze()

for site in mol_complex.interaction_sets.values():
    for interaction in site.all_itypes:
        # Objects are interaction-class-specific. Inspect the documented fields
        # for HydrophobicContact, HydrogenBond, PiStacking, SaltBridge, etc.;
        # there is no universal `.type` or `.ligatom.coords` interface.
        interaction_class = type(interaction).__name__
        print(interaction_class, interaction)
```

PLIP exposes typed interaction records with class-specific ligand/protein atoms and coordinates. Map those records to pharmacophore features explicitly and retain the interaction class and source atom identifiers.

## Pharmacophore Search (Pharmit / Pharmer)

For library screening, configure feature types, centers, radii, and optional shape constraints in Pharmit, or use a Pharmer database and query produced in the format required by the installed release. Do not pass LigandScout PML or a vendor `.ph4` file directly unless the selected interface documents that import path. Pharmit reported searching millions of conformers in seconds to minutes; actual runtime depends on query selectivity, database size, and deployment (Sunseri & Koes 2016).

## Pharmacophore Quality Validation

Evaluate a pharmacophore by:

1. **Retrospective enrichment**: a stated metric on target-relevant actives and inactives/decoys. DUD-E can provide a benchmark with known decoy-construction biases; COCONUT is a natural-products collection, not a target-specific active/decoy benchmark.
2. **Geometric tightness**: feature distance variance across actives
3. **Selectivity**: false positives in inactive set should be low
4. **Specific consistency**: pharmacophore matches each active's bioactive conformer

```python
def pharmacophore_enrichment(query_pharmacophore, actives, inactives,
                             matches_pharmacophore):
    """Return active/inactive match-rate enrichment for a supplied matcher."""
    if not actives or not inactives:
        raise ValueError('actives and inactives must both be non-empty')
    n_active_match = sum(
        bool(matches_pharmacophore(mol, query_pharmacophore))
        for mol in actives)
    n_inactive_match = sum(
        bool(matches_pharmacophore(mol, query_pharmacophore))
        for mol in inactives)
    active_rate = n_active_match / len(actives)
    inactive_rate = n_inactive_match / len(inactives)
    return float('inf') if inactive_rate == 0 else active_rate / inactive_rate
```

For this repository, enrichment >=5x may be used as a starting triage heuristic only after the active/decoy construction and matching policy are documented. Report the full metric and uncertainty, and calibrate the acceptance threshold on the project dataset.

## Pocket-Conditioned Pharmacophore Generation (PharmacoForge)

PharmacoForge (Flynn et al. 2025) applies a diffusion model to a protein pocket and generates candidate 3D pharmacophores. It does **not** directly generate molecular structures from an input pharmacophore. The validated workflow is:

1. Prepare the protein pocket in the representation required by the published PharmacoForge release.
2. Sample and rank pocket-conditioned pharmacophores.
3. Convert a selected pharmacophore into the query representation used by the search engine.
4. Retrieve matching, purchasable compounds and evaluate them with docking, strain, and physical-validity checks.

The paper compares pharmacophore and downstream retrieval performance with other pocket-based approaches; it does not support a drug-likeness or novelty comparison with REINVENT.

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

### PLIP -- water bridge absent from output

**Trigger:** Bridging water between ligand donor and protein acceptor.

**Mechanism:** PLIP can report water bridges, but the required crystallographic water must be present in the input and satisfy its geometric criteria.

**Symptom:** Pharmacophore missing critical H-bond feature.

**Fix:** Retain relevant crystallographic waters, inspect PLIP water-bridge output, and review borderline geometry manually.

## Reconciliation: Ligand-Based vs Receptor-Based

| Aspect | Ligand-based | Receptor-based |
|--------|--------------|----------------|
| Data needed | Multiple actives with defensible conformers/alignment | A defined pocket, optionally with a co-crystal ligand |
| Main bias | Known active chemotypes, conformer choice, and alignment | Pocket structure, protonation, retained waters, and interaction-detection/modeling rules |
| Hit-set behavior | Depends on feature abstraction and tolerances | Depends on selected pocket interactions, excluded volumes, and tolerances |
| Confidence evidence | Retrospective recovery across held-out actives/inactives | Recovery of known interaction geometry and retrospective or prospective validation |

Choose between ligand- and receptor-based models using the available structural/activity evidence and target-relevant validation. Neither approach is universally more reliable, diverse, or suitable for scaffold hopping.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `Pharm3D.EmbedPharmacophore` fails | Bounds matrix infeasible | Review/loosen justified bounds and, when more attempts are warranted, increase the documented `count` argument; inspect `n_failed` |
| Pharmacophore matches everything | Too few features | Add features; tighten tolerances |
| Pharmacophore matches nothing | Too many features or tight bounds | Reduce feature count; loosen tolerances |
| BaseFeatures.fdef not found | RDKit installation issue | Check `from rdkit.RDPaths import RDDataDir` |
| Pharmacophore-conformer mismatch | Wrong conformer used | Use bioactive conformer from crystal |
| Pharmit search timeout | Library too large | Pre-filter by 2D fingerprint Tanimoto |
| apo2ph4 PML has no useful model | No robust pocket hot spots at selected settings | Recheck pocket definition and documented thresholds; inspect alternative models |

## References

- Wolber & Langer, *J. Chem. Inf. Model.* 45:160-169 (2005) -- LigandScout pharmacophores. https://doi.org/10.1021/ci049885e
- Heider et al., *J. Chem. Inf. Model.* 63:101-110 (2023; published online 2022) -- apo2ph4. https://doi.org/10.1021/acs.jcim.2c00814
- Flynn EL, Shah R, Dunn I, Aggarwal R, Koes DR, *Front. Bioinform.* 5:1628800 (2025) -- PharmacoForge. https://doi.org/10.3389/fbinf.2025.1628800
- RDKit, `Chem.Pharm3D` API documentation. https://www.rdkit.org/docs/source/rdkit.Chem.Pharm3D.html
- RDKit, `EmbedPharmacophore` API documentation -- embedding molecules against an existing pharmacophore. https://www.rdkit.org/docs/source/rdkit.Chem.Pharm3D.EmbedLib.html#rdkit.Chem.Pharm3D.EmbedLib.EmbedPharmacophore
- COCONUT, official resource -- open natural-products collection. https://coconut.naturalproducts.net/
- Adasme et al., *Nucleic Acids Res.* 49:W530-W534 (2021) -- PLIP interaction profiler. https://doi.org/10.1093/nar/gkab294
- Sunseri & Koes, *Nucleic Acids Res.* 44:W442-W448 (2016) -- Pharmit interactive search. https://doi.org/10.1093/nar/gkw287

## Related Skills

- chemoinformatics/molecular-io - Parse molecules
- chemoinformatics/conformer-generation - Generate 3D for pharmacophore
- chemoinformatics/shape-similarity - 3D shape adjacent to pharmacophore
- chemoinformatics/virtual-screening - Pharmacophore as docking pre-filter
- chemoinformatics/scaffold-analysis - 2D scaffold-hopping context
- chemoinformatics/generative-design - Generate or optimize molecules after pharmacophore-based retrieval
- structural-biology/structure-io - PDB handling
