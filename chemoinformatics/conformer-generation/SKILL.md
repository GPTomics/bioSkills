---
name: bio-conformer-generation
description: Generates 3D conformer ensembles using RDKit ETKDGv3 with knowledge-enhanced distance geometry, MMFF94/UFF force-field optimization, CREST + GFN2-xTB semi-empirical refinement, and macrocycle-aware torsion preferences. Provides explicit decision rules for single vs ensemble conformer use, RMSD pruning, energy windows, conformer count, and force-field choice. Use when preparing 3D ligands for docking, generating descriptor input for 3D QSAR, or sampling macrocycle/peptide conformational ensembles.
tool_type: mixed
primary_tool: RDKit
---

## Version Compatibility

Reference examples tested with: RDKit 2024.09+, xtb 6.7+, CREST 3.0+, OpenMM 8.1+ for follow-up MD.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `xtb --version`; `crest --version`

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Conformer Generation

Generate 3D conformer ensembles for molecules from 2D structures. The choice of method depends on molecule size, flexibility, and downstream use: ETKDG / ETKDGv3 (Riniker & Landrum 2015 J Chem Inf Model 55:2562-2574; knowledge-enhanced distance geometry) is the modern default for drug-like molecules, MMFF94/UFF for fast energy minimization, CREST + GFN2-xTB for high-accuracy semi-empirical sampling of macrocycles and peptides. A single conformer is rarely sufficient: descriptor variance across the ensemble can exceed the descriptor signal, and docking pose accuracy degrades if the starting conformer is non-bioactive.

For docking pose validation, see `chemoinformatics/pose-validation`. For free-energy methods (which require ensemble sampling), see `chemoinformatics/free-energy-calculations`.

## Conformer Method Taxonomy

| Method | Cost / mol | Quality | Use case | Fails when |
|--------|-----------|---------|----------|------------|
| ETKDGv3 + MMFF94 | <1s | Good for drug-like | Default; docking input; descriptors | Macrocycles, peptides, transition metals |
| ETKDGv3 + UFF | <1s | Lower-quality MMFF94 alternative | Fallback when MMFF94 fails to parameterize | Same as MMFF94 |
| Omega (OpenEye) | 1s | Industry-standard commercial | Commercial pipelines | License cost |
| Confab (Open Babel) | 5s | Systematic torsion search | Patent expiration | Quality limited |
| RDKit ETKDGv3 + macrocycle preferences | 10-60s | Drug-like macrocycles | Macrocyclic peptides | Still limited; CREST better |
| CREST + GFN2-xTB | minutes | High-accuracy semi-empirical | Macrocycles, peptides, conformer ensembles for QSAR | Computationally expensive; metal centers |
| CREST + GFN-FF | seconds | GFN2 quality at FF speed | Quick screening | Limited element coverage |
| GeoMol (Ganea 2021) | <0.1s GPU | ML-fast, ETKDGv3-quality | Large library 3D conformers | ML training distribution |
| TorsionNet (Gogineni 2020) | <0.1s GPU | ML-fast | Drug-like | ML training distribution |
| MD sampling (OpenMM) | hours | High-quality dynamic | Free energy, induced fit | Computational cost |

**Decision:** For drug-like molecules (<500 Da, <8 rotatable bonds), **ETKDGv3 + MMFF94** with 20-100 conformers is the modern default. For macrocycles, peptides, or molecules with >12 rotatable bonds, **CREST + GFN2-xTB** captures the conformational diversity. For ML-scale (>1M molecules), **GeoMol** trades accuracy for speed.

## Decision Tree by Scenario

| Scenario | Method | Conformer count | Energy window |
|----------|--------|-----------------|----------------|
| Single docking pose (initial 3D) | ETKDGv3 + MMFF94 | 1 | n/a |
| Multi-conformer docking | ETKDGv3 + MMFF94 | 10-50 | 10 kcal/mol |
| 3D QSAR descriptor input | ETKDGv3 + MMFF94 | 50-200 | 5 kcal/mol |
| Pharmacophore search | ETKDGv3 + MMFF94 | 100-500 | 5 kcal/mol |
| Macrocycle / peptide | CREST + GFN2-xTB | 50-200 (auto from CREST) | 5-8 kcal/mol |
| FEP input | CREST + GFN2-xTB then MD relax | 1-3 representative | 3 kcal/mol |
| Bioactive conformer search | ETKDGv3 + MMFF94 then dock with rescore | 100-500 | 10 kcal/mol |
| Shape similarity / ROCS | ETKDGv3 + MMFF94 | 50-200 | 10 kcal/mol |
| Conformer-dependent descriptors | ETKDGv3 ensemble + Boltzmann avg | 20-100 | 5 kcal/mol |

## ETKDGv3 (Modern Default)

ETKDGv3 (Riniker & Landrum 2015) incorporates experimental torsion preferences into distance geometry: starts from random embeddings, refines by satisfying experimentally-derived bond, angle, and torsion preferences.

**Goal:** Generate an ensemble of 3D conformers from a SMILES with the modern default embedding algorithm.

**Approach:** Add explicit hydrogens, configure ETKDGv3 params (random seed, max attempts, random coords), and embed multiple conformers via `EmbedMultipleConfs`.

```python
from rdkit import Chem
from rdkit.Chem import AllChem

def gen_conformers(smiles, n_conf=20, seed=42):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.useRandomCoords = True
    params.maxAttempts = 1000
    ids = AllChem.EmbedMultipleConfs(mol, numConfs=n_conf, params=params)
    return mol, list(ids)
```

`useRandomCoords=True` improves convergence for macrocycles and heavily-rotated molecules. `maxAttempts=1000` handles difficult embeddings.

## Force-Field Optimization

After embedding, minimize each conformer to a local minimum.

**Goal:** Reduce strain in each embedded conformer to a stable local minimum and record the resulting energies.

**Approach:** Build MMFF94s force-field parameters, minimize each conformer in place, and collect energies; fall back to UFF when MMFF94 cannot parameterize the molecule.

```python
def optimize_conformers(mol, conf_ids, force_field='mmff94'):
    energies = []
    if force_field == 'mmff94':
        mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')
        for cid in conf_ids:
            ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=cid)
            ff.Minimize()
            energies.append(ff.CalcEnergy())
    else:  # UFF fallback
        for cid in conf_ids:
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=cid)
            ff.Minimize()
            energies.append(ff.CalcEnergy())
    return energies
```

**MMFF94 vs MMFF94s:** MMFF94s is the "standard" set with simpler aromatic nitrogen handling; preferred for most drug-like.

**UFF (Universal Force Field):** Lower quality but handles any element including transition metals. Use as fallback when MMFF94 cannot parameterize (uncommon elements, charged species).

## RMSD Pruning

Remove near-duplicate conformers within a chosen RMSD cutoff to keep the ensemble diverse:

```python
import numpy as np

def prune_conformers_rmsd(mol, conf_ids, rmsd_cutoff=0.5):
    n = len(conf_ids)
    keep = []
    for i, cid in enumerate(conf_ids):
        is_unique = True
        for kept_cid in keep:
            rmsd = AllChem.GetBestRMS(mol, mol, cid, kept_cid)
            if rmsd < rmsd_cutoff:
                is_unique = False
                break
        if is_unique:
            keep.append(cid)
    return keep
```

**Typical RMSD cutoff (Source / Rationale):**

| Cutoff | Use case | Source |
|--------|----------|--------|
| 0.5 Å | Drug-like ensemble for descriptors / docking | Empirical: below this conformers represent same minimum (Hawkins 2007) |
| 1.0 Å | Drug-like ensemble for pharmacophore | Standard ROCS / pharmacophore practice |
| 1.5-2.0 Å | Macrocycles / peptides | Higher conformational freedom; Tan 2018 macrocycle benchmarks |
| 2.0+ Å | Cluster-centroid representative ensembles | Coarse representative sampling |

## Energy Window Filtering

Remove conformers above an energy cutoff (high-energy conformers are unlikely to be bioactive):

```python
def filter_by_energy(mol, conf_ids, energies, window_kcal=10.0):
    min_e = min(energies)
    keep = []
    for cid, e in zip(conf_ids, energies):
        if e - min_e <= window_kcal:
            keep.append(cid)
    return keep
```

**Window choice:**
- 3 kcal/mol: very strict, only near-global-min conformers (FEP, MD setup)
- 5 kcal/mol: typical for 3D QSAR, pharmacophore
- 10 kcal/mol: typical for docking input (bioactive conformer may be higher)
- 25 kcal/mol: macrocycles, no filter (bioactive conformer can be high-energy when bound)

## Macrocycle Handling

Macrocycles (>=12 atom rings) have distinct conformational issues: ETKDGv3 default knowledge base under-samples macrocycle torsions. Use macrocycle-specific torsion preferences:

```python
from rdkit.Chem import AllChem

def macrocycle_conformers(smiles, n_conf=200, seed=42):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.useRandomCoords = True
    params.useMacrocycleTorsions = True
    params.useSmallRingTorsions = True
    params.maxAttempts = 5000
    ids = AllChem.EmbedMultipleConfs(mol, numConfs=n_conf, params=params)
    return mol, list(ids)
```

For pharmaceutical macrocycles (cyclosporine, paclitaxel, large peptides), CREST + GFN2-xTB is the gold standard.

## CREST + GFN2-xTB for High-Quality Sampling

CREST (Grimme 2024) performs iterative meta-dynamics + GFN2-xTB optimization for conformer sampling.

**Goal:** Sample high-quality conformer ensembles for macrocycles, peptides, or molecules where ETKDGv3 + MMFF94 is inadequate.

**Approach:** Start from an RDKit-generated MMFF94-relaxed conformer, write to XYZ, and run CREST with GFN2-xTB driver to perform iterative meta-dynamics + reoptimization.

```bash
xtb mol.xyz --opt extreme
crest opt.xyz --gfn2 --T 12 -ewin 6
```

**`--gfn2`**: use GFN2-xTB (most accurate of GFN family for drug-like molecules).
**`--gfn-ff`**: use GFN-FF (faster, less accurate).
**`-ewin 6`**: 6 kcal/mol energy window above global min.
**`-T 12`**: use 12 CPU threads.

Output: `crest_conformers.xyz` with sampled ensemble.

**Workflow:** Start from RDKit ETKDGv3 + MMFF94 (cheap initial structure) -> save as XYZ -> CREST refinement.

```python
from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess

def crest_workflow(smiles, out_dir='crest_out'):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol)

    xyz = Chem.MolToXYZBlock(mol)
    with open(f'{out_dir}/input.xyz', 'w') as f:
        f.write(xyz)
    subprocess.run(['crest', f'{out_dir}/input.xyz', '--gfn2', '-T', '12'],
                   cwd=out_dir, check=True)
    return f'{out_dir}/crest_conformers.xyz'
```

## Boltzmann Averaging of Properties

For ensemble descriptors (3D shape, dipole moment, polar surface area in 3D), Boltzmann-weight by energy:

```python
import numpy as np

def boltzmann_weights(energies, T=300.0):
    energies = np.array(energies)
    kt = 0.001987 * T  # kcal/mol at 300K
    rel = energies - energies.min()
    w = np.exp(-rel / kt)
    return w / w.sum()

def boltzmann_average(values, energies, T=300.0):
    w = boltzmann_weights(energies, T)
    return float(np.sum(np.array(values) * w))
```

For Boltzmann averaging, energies should be MMFF94 or higher quality. UFF energies are unreliable for Boltzmann weighting.

## ML-Based Conformer Generation (GeoMol, TorsionNet)

For very large libraries (>1M compounds), classical methods become bottlenecks. ML-based methods generate conformers in <0.1s/mol on GPU:

```python
# Pseudo-code for GeoMol-style ML conformer generation
# (Requires pre-trained model + dependencies)
# from geomol import generate_conformers
# conformers = generate_conformers(smiles, n_conformers=10)
```

**Trade-off:** ML methods (GeoMol, TorsionNet) match ETKDGv3 quality on drug-like molecules but extrapolate poorly outside training distribution (macrocycles, organometallics).

## Per-Tool Failure Modes

### ETKDGv3 -- failed embedding

**Trigger:** Macrocycle, highly constrained polycyclic, or sterically crowded molecule.

**Mechanism:** Distance geometry cannot find consistent 3D structure within max attempts.

**Symptom:** `EmbedMolecule` returns -1; `EmbedMultipleConfs` returns empty list.

**Fix:** Set `useRandomCoords=True`, increase `maxAttempts` to 5000+; for macrocycles, set `useMacrocycleTorsions=True`. As fallback, use CREST.

### MMFF94 -- parameter missing

**Trigger:** Molecule contains element not parameterized (transition metals, certain S+ species).

**Mechanism:** MMFF94 only covers H, C, N, O, F, Si, P, S, Cl, Br, I + select cations.

**Symptom:** `MMFFGetMoleculeProperties` returns None; optimization silently no-ops.

**Fix:** Fall back to UFF; or for metals, use GFN2-xTB.

### Conformer ensemble too small

**Trigger:** `n_conf=10` for a flexible molecule (>5 rotatable bonds).

**Mechanism:** 10 conformers insufficient to sample conformational space; many minima missed.

**Symptom:** RMSD distribution narrow; descriptor variance underestimated.

**Fix:** Use n_conf = max(10, 5 * NumRotatableBonds + 10) heuristic (Hawkins 2017).

### Single-conformer 3D descriptor

**Trigger:** Calculating 3D descriptors from a single conformer.

**Mechanism:** 3D descriptor variance across conformers can be 50%+ of mean.

**Symptom:** Same molecule produces different 3D descriptors on rerun.

**Fix:** Always compute descriptor over ensemble; report mean ± std, or Boltzmann-weighted mean.

### CREST -- timeout on flexible molecule

**Trigger:** Cyclosporin or large peptide.

**Mechanism:** CREST metadynamics scales poorly with rotational complexity.

**Symptom:** Hours of CPU time per molecule; incomplete sampling.

**Fix:** Use `--gfn-ff` for faster initial sampling; reduce metadynamics time `--mdtime 5` or skip metadyn with `--noopt`.

### GFN2-xTB conformer reordering

**Trigger:** Comparing conformer energies between GFN2-xTB and DFT.

**Mechanism:** GFN2-xTB is parameterized for energies; relative conformer ordering can differ from DFT by 1-2 kcal/mol.

**Symptom:** "Wrong" conformer reported as global minimum vs DFT reference.

**Fix:** For high-stakes work, re-rank top GFN2-xTB conformers with DFT single-points (e.g., r2SCAN-3c).

## Reconciliation: ETKDGv3 vs CREST

| Use case | ETKDGv3 | CREST |
|----------|---------|-------|
| Drug-like, <500 Da, <8 RotBonds | Sufficient | Overkill |
| 8-12 RotBonds | OK with n_conf>=100 | Better at expense of cost |
| Macrocycle, peptide, >12 RotBonds | Inadequate | Required |
| Boltzmann-weighted descriptors | OK but energies less accurate | Better |
| FEP input | Possible | Preferred (after MMFF cleanup) |

For ETKDGv3 ensembles, run CREST on a subset for benchmarking; if RMSD < 1A across methods, ETKDGv3 is adequate.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `EmbedMolecule` returns -1 | Embed failed | Set `useRandomCoords=True`; raise `maxAttempts` |
| MMFFOptimize no-op | MMFF parameters missing | Use UFF fallback |
| All conformers identical | Stiff molecule | OK; molecule is rigid |
| Conformers physically wrong | Stereochemistry lost | Re-add explicit stereo before embedding |
| 3D descriptors differ per run | Random seed not set | `params.randomSeed = 42` |
| CREST out-of-memory | Too many conformers in search | Reduce `--T` threads; raise `--ewin` window |
| Macrocycle ring inverted | Default torsion preferences wrong | Set `useMacrocycleTorsions=True` |
| AddHs not called | Implicit H not embedded | `mol = Chem.AddHs(mol)` before EmbedMolecule |

## References

- Hawkins et al., *J. Chem. Inf. Model.* 50:572 (2010) -- OMEGA conformer sampling (commercial OpenEye tool; predecessor to ETKDG).
- Riniker & Landrum, *J. Chem. Inf. Model.* 55:2562-2574 (2015) -- ETKDG / ETKDGv3 with torsion preferences (the standard RDKit knowledge-distance-geometry embedder).
- Halgren, *J. Comput. Chem.* 17:490 (1996) -- MMFF94 force field.
- Rappe et al., *J. Am. Chem. Soc.* 114:10024 (1992) -- UFF.
- Pracht, Bohle, Grimme, *J. Chem. Phys.* 160:114110 (2024) -- CREST 3.0.
- Bannwarth et al., *J. Chem. Theory Comput.* 15:1652 (2019) -- GFN2-xTB.
- Ganea et al., *NeurIPS* (2021) -- GeoMol ML conformer generation.
- Hawkins, *J. Chem. Inf. Model.* 57:1747 (2017) -- conformer count heuristics.

## Related Skills

- chemoinformatics/molecular-io - Parse molecules
- chemoinformatics/molecular-standardization - Standardize before embedding
- chemoinformatics/molecular-descriptors - 3D descriptors from ensembles
- chemoinformatics/shape-similarity - Multi-conformer 3D shape matching
- chemoinformatics/virtual-screening - Generate 3D ligands for docking
- chemoinformatics/free-energy-calculations - Sample conformers for MD setup
- chemoinformatics/pharmacophore-modeling - 3D pharmacophore from ensembles
