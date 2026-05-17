---
name: bio-shape-similarity
description: Performs 3D shape-based similarity searching using ROCS (OpenEye), USRCAT (ultra-fast), Open3DAlign (RDKit), ESPSim (electrostatic), and ShaEP with explicit handling of Tanimoto-Combo (shape + color), shape vs ECFP4 complementarity, conformer-ensemble searching, alignment optimization, and scaffold hopping. Use when searching for shape-mimicking compounds with different scaffolds, identifying bioisosteric replacements, prospective scaffold hopping, or expanding hit series beyond 2D similarity.
tool_type: python
primary_tool: RDKit
---

## Version Compatibility

Reference examples tested with: RDKit 2024.09+ (Open3DAlign), USRCAT 1.2+, ShaEP 1.7+, ROCS (OpenEye, commercial).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Shape Similarity

Search for compounds with similar 3D shape (and optionally chemical features) to a query molecule. Shape-based screening complements 2D fingerprint search: it can find scaffold-hopped compounds that ECFP4 misses (different scaffolds with similar shape). ROCS (OpenEye) is the industry-standard commercial tool; Open3DAlign (RDKit), USRCAT (Schreyer & Blundell 2012), and ShaEP are open-source alternatives. Modern best practice combines shape with color (chemical-feature similarity) via Tanimoto-Combo: matches share both shape and pharmacophore feature distribution.

For 2D fingerprint similarity, see `chemoinformatics/similarity-searching`. For pharmacophore search (discrete feature constraints), see `chemoinformatics/pharmacophore-modeling`. For 3D conformer generation, see `chemoinformatics/conformer-generation`.

## Shape Method Taxonomy

| Tool | Speed | Approach | Open-source | Fails when |
|------|-------|----------|-------------|------------|
| ROCS (OpenEye) | 1k mols/sec on GPU (FastROCS) | Gaussian shape + color | No | License cost |
| ROCS X (Sept 2025) | Multi-billion library, GPU | ML-enhanced shape | No | Limited release |
| USRCAT | 100k mols/sec | Ultrafast moment-based + atom types | Yes | Coarse approximation |
| Open3DAlign (RDKit) | 100 mols/sec | Iterative volume overlap | Yes | Optimization slow |
| ShaEP | 10 mols/sec | Field-based (shape + ESP) | Yes | Less standard |
| ESPSim | similar to ShaEP | Electrostatic + shape | Yes | Limited public benchmarks |
| Phase-Shape (Schrödinger) | commercial | Shape + pharmacophore | No | Commercial |
| USR (original) | 100k mols/sec | Moment-based only | Yes | No atom type info |

**Decision:** For commercial pipelines, **ROCS** is the gold standard. For open-source, **Open3DAlign** is the most accurate; **USRCAT** is the fastest for ultralarge libraries.

## Decision Tree by Scenario

| Scenario | Method | Notes |
|----------|--------|-------|
| Lead-like library, search top 100k | USRCAT pre-filter + Open3DAlign rescore | Hybrid speed/accuracy |
| Production VS for scaffold hop | ROCS + color (commercial) | Industry standard |
| Scaffold hopping prospective | Open3DAlign with conformer ensemble | Shape + flexibility |
| Bioisostere replacement | ROCS color with neutral scoring | Pharmacophore-equivalent matches |
| Patent space carve-out | Shape constraint + 2D dissimilarity | Combine shape + dissimilar scaffold |
| Library diversity assessment | USRCAT k-nearest neighbor | Fast |
| Crystal-bound conformer template | Open3DAlign starting from co-crystal pose | Bioactive shape |
| Cross-target screening | Shape + pharmacophore feature | Combined screen |

## Tanimoto-Combo Scoring (ROCS Standard)

Tanimoto-Combo = (Tanimoto_shape + Tanimoto_color) / 2

- Tanimoto_shape: volume overlap normalized
- Tanimoto_color: pharmacophore feature overlap

| Range | Interpretation |
|-------|----------------|
| > 1.0 | Very similar shape + color (rare; top hits) |
| 0.7-1.0 | Strong hit; likely binding mode similarity |
| 0.5-0.7 | Moderate; further validation needed |
| 0.3-0.5 | Weak; many false positives |
| < 0.3 | Background |

In ROCS production, hits with TanimotoCombo > 0.7 are typically followed up.

## USRCAT (Ultra-Fast Shape Recognition + Atom Types)

USRCAT (Schreyer & Blundell 2012) extends Ultrafast Shape Recognition (USR) with atom-type information. Each molecule is represented as a 60-dimensional moment vector (12 moments × 5 atom types).

**Goal:** Encode a molecule into the 60-D USRCAT moment vector and score similarity against another molecule for alignment-free shape search.

**Approach:** Parse the SMILES, add hydrogens, generate one 3D conformer with ETKDGv3, compute USRCAT descriptors, and apply the inverse-mean-absolute-difference similarity to a second descriptor vector.

```python
from usrcat import compute_usrcat_descriptors, compute_similarity

mol = Chem.MolFromSmiles('CCO')
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())

descriptors = compute_usrcat_descriptors(mol)
# Returns numpy array of 60 floats: 12 USR moments x 5 atom types
# (hydrophobic, aromatic, acceptor, donor, anion/cation)
# Similarity between two descriptor vectors: 1 / (1 + mean_abs_difference)

similarity = compute_similarity(desc1, desc2)  # 0-1, higher = more similar
```

**Speed:** O(N) descriptor calculation (no alignment); O(1) similarity comparison. Suitable for >10M compound libraries.

**Limit:** USRCAT is a coarse approximation. Predictive for analog identification; less precise for scaffold hopping.

## Open3DAlign (RDKit)

Open3DAlign performs iterative alignment to maximize volume overlap:

**Goal:** Align a target molecule onto a query in 3D and score volume overlap with Open3DAlign.

**Approach:** Build 3D structures for query and target (parse SMILES, add hydrogens, ETKDG embed), run `GetO3A` to find the best alignment, then call `Align()` for the in-place RMSD and `Score()` for the overlap score.

```python
from rdkit.Chem import rdMolAlign

query = Chem.MolFromSmiles('CCC(=O)Nc1ccccc1')
query = Chem.AddHs(query)
AllChem.EmbedMolecule(query, AllChem.ETKDGv3())

target = Chem.MolFromSmiles('CCC(=O)Nc1ccc(F)cc1')
target = Chem.AddHs(target)
AllChem.EmbedMolecule(target, AllChem.ETKDGv3())

O3A = rdMolAlign.GetO3A(target, query)
rmsd = O3A.Align()  # aligns target to query in place
score = O3A.Score()
```

`GetO3A` finds best alignment between conformers; `Align()` aligns and returns RMSD; `Score()` returns Open3DAlign score (similar to TanimotoCombo).

**Open3DAlign vs ROCS:** Open3DAlign is open-source and competitive on small benchmarks; slower than ROCS at scale.

## Conformer-Ensemble Shape Searching

For each library molecule, generate ensemble of conformers; pick best-shape conformer:

**Goal:** Run shape-similarity search over a conformer ensemble per library molecule so bound-conformer-like shapes are recovered.

**Approach:** For each library molecule, add hydrogens, embed n_conf conformers with ETKDGv3, MMFF-optimize, score each conformer against the query with Open3DAlign, and keep the best score per molecule.

```python
def shape_search_ensemble(query_mol, library_mols, n_conf=20):
    hits = []
    for target in library_mols:
        target = Chem.AddHs(target)
        AllChem.EmbedMultipleConfs(target, numConfs=n_conf,
                                    params=AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMoleculeConfs(target)

        scores = []
        for c in range(target.GetNumConformers()):
            O3A = rdMolAlign.GetO3A(target, query_mol, prbCid=c)
            scores.append(O3A.Score())
        hits.append((target, max(scores)))
    return sorted(hits, key=lambda x: x[1], reverse=True)
```

**Critical:** Single-conformer shape search misses ~30% of true hits because the wrong conformer is sampled. Always use ensemble.

## ESP Similarity (Electrostatic)

ShaEP and ESPSim extend shape with electrostatic surface potential overlap. For ESP-relevant pharmacophores (binding pockets with strong electrostatics):

```bash
shaep --query query.mol2 --target target.mol2 --output match.sdf \
      --esp-weight 0.5
```

ESP scoring catches electrostatic-equivalent bioisosteres that pure shape misses (carboxylate vs tetrazole same charge).

## Shape vs ECFP4 Complementarity

| Shape Tanimoto | ECFP4 Tanimoto | Interpretation |
|----------------|----------------|----------------|
| > 0.7 | > 0.7 | Same chemotype, same shape (close analog) |
| > 0.7 | < 0.5 | Scaffold-hop! Different chemotype, similar shape |
| < 0.5 | > 0.7 | Same chemotype, different shape (flexible) |
| < 0.5 | < 0.5 | Unrelated |

The shape >> ECFP4 quadrant is the scaffold-hopping gold:

**Goal:** Identify scaffold-hop candidates that are 3D-shape-similar but 2D-chemotype-dissimilar to the query.

**Approach:** Run the conformer-ensemble shape search, keep hits above a shape Tanimoto cutoff, then retain only those whose ECFP4 Tanimoto to the query is below an ECFP4 dissimilarity cutoff.

```python
def scaffold_hop_candidates(query_mol, library, shape_threshold=0.7,
                            ecfp_threshold=0.5):
    shape_hits = shape_search_ensemble(query_mol, library)
    candidates = []
    for target, shape_score in shape_hits:
        if shape_score >= shape_threshold:
            ecfp_sim = ecfp_tanimoto(query_mol, target)
            if ecfp_sim < ecfp_threshold:
                candidates.append((target, shape_score, ecfp_sim))
    return candidates
```

## Per-Tool Failure Modes

### USRCAT -- false positive on small molecules

**Trigger:** Library has many fragment-sized compounds.

**Mechanism:** USRCAT moments dominated by overall shape; small molecules look "similar" if shape resemble.

**Symptom:** Many fragment hits; not pharmacophore-relevant.

**Fix:** Filter by MW (>= 200); use Open3DAlign for rescoring.

### Open3DAlign -- slow on large library

**Trigger:** Million-compound library, full alignment.

**Mechanism:** Open3DAlign is iterative; O(N) per molecule.

**Symptom:** Hours of compute.

**Fix:** Pre-filter with USRCAT (fast), Open3DAlign on top 10k.

### Shape only -- wrong stereochemistry match

**Trigger:** Mirror-image of correct binder.

**Mechanism:** Pure shape Tanimoto symmetric under chirality inversion.

**Symptom:** Enantiomer of inactive scores as hit.

**Fix:** Validate hits by 3D pose; check stereochemistry.

### ROCS color -- bioisostere missed

**Trigger:** -COOH replaced by -SO3H or tetrazole.

**Mechanism:** Default color types may not equate these bioisosteres.

**Symptom:** Known bioisostere doesn't score high.

**Fix:** Use color-only scoring without size penalty; or pharmacophore-feature-equivalence.

### Conformer not bioactive

**Trigger:** Library compound generated conformer is not the bound conformation.

**Mechanism:** ETKDGv3 generates plausible conformers; bound conformer may be higher energy.

**Symptom:** Known active doesn't shape-match query.

**Fix:** Use larger conformer ensemble; weight by Boltzmann; or use CREST + GFN2-xTB for high-quality sampling.

### Field-based methods slower

**Trigger:** ShaEP or ESPSim on production library.

**Mechanism:** Field-based methods compute Gaussian fields per molecule.

**Symptom:** 10-100x slower than ROCS shape-only.

**Fix:** Use as second-stage rescore; not primary screen.

## Reconciliation: Shape vs Pharmacophore

| Aspect | Shape | Pharmacophore |
|--------|-------|----------------|
| Representation | Volume distribution | Discrete features in space |
| Captures | Overall bulk | Interaction-relevant features |
| Speed | Fast (USRCAT) to medium (Open3DAlign) | Fast |
| Specificity | Lower | Higher |
| False positive rate | Medium | Lower |
| Best for | Scaffold hopping initial | Scaffold hopping refinement |

Use shape for *broad search* (high recall, moderate precision); use pharmacophore for *refinement* (lower recall, high precision).

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Open3DAlign returns 0 | No reasonable alignment found | Use random starting rotation; increase attempts |
| USRCAT vector all zeros | Mol has no 3D coords | Generate conformer first |
| Shape Tanimoto > 1 | Implementation bug or unnormalized | Check formula; ROCS reports unnormalized possible |
| ROCS very slow | Sequential processing | Use parallel batching |
| Shape match but no docking pose | Wrong binding pose | Use docking on top shape hits, not shape alone |
| Missing co-crystal template | Apo or AlphaFold-only structure | Use ligand-based pharmacophore + shape |
| ShaEP returns no hits | Strict tolerance | Loosen overlap thresholds |

## References

- Hawkins et al., *J. Med. Chem.* 50:74 (2007) -- ROCS algorithm.
- Schreyer & Blundell, *J. Cheminformatics* 4:27 (2012) -- USRCAT.
- Vainio & Johnson, *J. Chem. Inf. Model.* 47:2462 (2007) -- ShaEP.
- Liu et al., *J. Chem. Inf. Model.* (2024) -- Open3DAlign improvements.
- Roy et al., *J. Med. Chem.* 65:11875 (2022) -- shape-based VS modern review.

## Related Skills

- chemoinformatics/molecular-io - Parse query and library
- chemoinformatics/conformer-generation - Generate 3D conformer ensembles
- chemoinformatics/similarity-searching - 2D similarity comparison
- chemoinformatics/pharmacophore-modeling - Pharmacophore alternative
- chemoinformatics/scaffold-analysis - 2D scaffold analysis
- chemoinformatics/virtual-screening - Shape as pre-filter to docking
