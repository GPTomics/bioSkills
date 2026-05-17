---
name: bio-pose-validation
description: Validates docked / generated protein-ligand poses using PoseBusters physical-validity tests, strain energy quantification, geometric checks (planarity, vdW overlap, bond/angle distortion), and pose-energy reasonableness. Filters AI-docking outputs (DiffDock, EquiBind, NeuralPLexer) where ~50% of poses fail physical-validity tests. Use when QC-ing docking results, comparing classical vs ML docking outputs, or filtering pose lists before SAR analysis.
tool_type: python
primary_tool: PoseBusters
---

## Version Compatibility

Reference examples tested with: PoseBusters 0.6+, RDKit 2024.09+, pandas 2.2+, posecheck 0.5+ (optional).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Pose Validation

Test docked or AI-generated protein-ligand poses for physical plausibility. PoseBusters (Buttenschoen 2024) is the modern gold standard: a suite of geometric, chemical, and energetic checks that flag implausible poses (planar aromatic rings now non-planar, vdW clashes, broken bonds, wrong chirality, unrealistic torsions). The PoseBusters benchmark showed that AI-based docking methods (DiffDock, EquiBind, TANKBind) produce ~50% physically-invalid poses despite reporting good RMSD; classical methods (Vina, GOLD) produce ~5-15% invalid. PB-valid status is therefore essential for downstream SAR, FEP setup, or generative model training.

For docking, see `chemoinformatics/virtual-screening`. For ML docking specifically, see `chemoinformatics/ml-docking-rescoring`.

## PoseBusters Test Suite

PoseBusters runs ~20 individual checks grouped into:

| Check group | What it tests | Threshold |
|-------------|---------------|-----------|
| Sanity | Ligand chemical sanity | RDKit sanitization passes |
| Bond lengths | Bond lengths within reference | < 2 std from RDKit defaults |
| Bond angles | Bond angles within reference | < 2 std from RDKit defaults |
| Internal steric | No intra-ligand vdW clash | vdW overlap < 1.0 Å |
| Aromatic ring planarity | Aromatic rings planar | < 0.25 Å RMS deviation |
| Double-bond stereo | Z/E preserved | Match input SMILES |
| Internal energy | Strain not absurd | UFF energy < 100 kcal/mol typical |
| Volume overlap | vdW overlap with protein | < 7.5% of ligand vdW volume |
| Distance to protein | Not floating in solvent | Closest contact < 5 Å |
| Chirality | R/S preserved from input | Match input SMILES |

A pose passing ALL tests is "PB-valid". Combined PB-valid + RMSD <= 2 Å is the modern criterion.

## When to Apply PoseBusters

| Workflow | PoseBusters use | Action |
|----------|-----------------|--------|
| Self-docking (validating method) | Required | Compare PB-valid + RMSD <= 2A |
| Cross-docking | Required | PB-valid + RMSD <= 2A; account for protein flexibility |
| Virtual screening top hits | Required | Filter to PB-valid before MM/GBSA / FEP |
| AI docking (DiffDock, etc.) | Mandatory | Often 50% fail; otherwise can't compare to classical |
| Generated structures (RFdiffusion) | Required | Generation often produces clashes |
| Boltz-2 / AlphaFold3 ligand poses | Recommended | Same family of issues; less frequent than DiffDock |
| Production FEP setup | Required | Strain in input pose breaks FEP convergence |

## PoseBusters Usage

```python
from posebusters import PoseBusters

bust = PoseBusters(config='redock')

results = bust.bust(
    mol_pred='predicted.sdf',
    mol_true='reference.sdf',
    mol_cond='receptor.pdb',
)
```

`config` options and included checks:

| Config | Includes | When to use |
|--------|----------|-------------|
| `redock` | All checks + RMSD vs reference + protein vdW overlap | Self-docking benchmarks, retrospective validation |
| `dock` | All checks except RMSD reference | Blind docking, prospective virtual screening |
| `mol` | Intra-ligand only (sanity, bonds, angles, rings, stereo, energy) | Conformer QC; no protein context |

Output DataFrame columns include: `mol_pred_loaded`, `sanitization`, `all_atoms_connected`, `bond_lengths`, `bond_angles`, `internal_steric_clash`, `aromatic_ring_flatness`, `double_bond_flatness`, `internal_energy`, `protein_flexibility`, `minimum_distance_to_protein`, `minimum_distance_to_organic_cofactors`, `minimum_distance_to_inorganic_cofactors`, `volume_overlap_with_protein`, `volume_overlap_with_organic_cofactors`, `volume_overlap_with_inorganic_cofactors`. All bool; True = pass.

Output: DataFrame with one column per check, one row per pose, boolean pass/fail.

## Python Library API

**Goal:** Programmatically validate a docked-pose SDF against a receptor PDB and produce a PB-valid filter.

**Approach:** Instantiate `PoseBusters(config='dock')`, call `bust()` on the SDF + PDB pair, and AND-aggregate all boolean check columns into a single `pb_valid` flag.

```python
from posebusters import PoseBusters
import pandas as pd

bust = PoseBusters(config='dock')

results = bust.bust(
    mol_pred='/path/to/docked_poses.sdf',
    mol_cond='/path/to/receptor.pdb',
)

results['pb_valid'] = results.iloc[:, 4:].all(axis=1)
valid = results[results['pb_valid']]
print(f'{len(valid)} / {len(results)} poses are PB-valid')
```

## Strain Energy Quantification

Beyond binary PB-valid, quantitative strain energy distinguishes "marginal" from "egregious" poses.

**Goal:** Quantify how far each docked pose is from its lowest-energy free conformer in MMFF94 energy units.

**Approach:** Generate a reference conformer ensemble (ETKDGv3 + MMFF94), take the global minimum energy as baseline, and report `docked_energy - min_ref_energy` as strain.

```python
from rdkit import Chem
from rdkit.Chem import AllChem

def ligand_strain(docked_sdf, n_ref=20):
    suppl = Chem.SDMolSupplier(docked_sdf, removeHs=False)
    strains = []
    for docked in suppl:
        if docked is None:
            continue

        smi = Chem.MolToSmiles(docked)
        ref = Chem.MolFromSmiles(smi)
        ref = Chem.AddHs(ref)
        AllChem.EmbedMultipleConfs(ref, numConfs=n_ref, params=AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMoleculeConfs(ref)

        ref_energies = []
        for c in range(ref.GetNumConformers()):
            ff = AllChem.MMFFGetMoleculeForceField(
                ref,
                AllChem.MMFFGetMoleculeProperties(ref),
                confId=c
            )
            ref_energies.append(ff.CalcEnergy())
        min_ref = min(ref_energies)

        docked_ff = AllChem.MMFFGetMoleculeForceField(
            docked,
            AllChem.MMFFGetMoleculeProperties(docked)
        )
        docked_e = docked_ff.CalcEnergy() if docked_ff else None

        strains.append({
            'min_ref_energy': min_ref,
            'docked_energy': docked_e,
            'strain': docked_e - min_ref if docked_e else None,
        })
    return strains
```

**Strain interpretation:**
- < 2 kcal/mol: ideal; close to relaxed
- 2-5 kcal/mol: realistic for bioactive conformer
- 5-10 kcal/mol: marginal; check for plausibility
- > 10 kcal/mol: implausible without exceptional binding

(Boström 1998 showed bioactive ligands have median strain ~2-4 kcal/mol.)

## vdW Overlap with Protein

PoseBusters limit: protein-ligand vdW overlap < 7.5% of ligand vdW volume. Computed via:

```python
def vdw_overlap_fraction(ligand_pdb, protein_pdb):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    # Use RDKit AllChem.UFFGetMoleculeForceField for vdW radii
    # Compute pairwise distances, sum overlap volume
    # (full implementation in PoseBusters; this is a sketch)
    return overlap_fraction
```

In practice, PoseBusters' `bust(...)` does this automatically.

## Aromatic Ring Planarity

```python
import numpy as np

def aromatic_planarity(mol):
    deviations = []
    for ring in mol.GetRingInfo().AtomRings():
        ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if not all(a.GetIsAromatic() for a in ring_atoms):
            continue
        coords = np.array([mol.GetConformer().GetAtomPosition(i)
                          for i in ring])
        centroid = coords.mean(axis=0)
        centered = coords - centroid
        _, s, vh = np.linalg.svd(centered)
        normal = vh[-1]
        deviation = np.abs(centered @ normal).max()
        deviations.append(deviation)
    return max(deviations) if deviations else 0
```

Aromatic ring deviation > 0.25 Å is implausible; flag.

## Per-Tool Failure Modes

### DiffDock-L -- chirality inversion

**Trigger:** Default DiffDock-L on chiral ligand.

**Mechanism:** Diffusion sampling can invert stereocenters; not constrained by training.

**Symptom:** PoseBusters chirality check fails; SMILES from pose != input SMILES.

**Fix:** Filter to PB-valid; rerun with constrained stereo if supported.

### EquiBind -- planar aromatic violation

**Trigger:** EquiBind single-shot pose prediction.

**Mechanism:** Equivariant NN does not preserve ring planarity.

**Symptom:** Aromatic rings buckled.

**Fix:** Post-relax with UFF/MMFF94 minimization before downstream use.

### TANKBind -- vdW clash

**Trigger:** TANKBind pose in tight pocket.

**Mechanism:** Distance prediction not constrained to vdW radii.

**Symptom:** PoseBusters vdW overlap fails.

**Fix:** Post-process with constrained energy minimization (UFF + frozen heavy atoms 10% relaxation).

### Boltz-1 -- bond length distortion

**Trigger:** Cofactor-bound complex.

**Mechanism:** Model trained on bond-length distribution; uncommon ligand pushes outside training distribution.

**Symptom:** Bond lengths >2 std from RDKit reference.

**Fix:** Post-relax with MMFF94 (constrained to maintain pose).

### High strain after Vina docking

**Trigger:** Highly constrained pocket; flexible ligand.

**Mechanism:** Vina exhaustiveness limits sampling; not enough rotation search.

**Symptom:** Strain energy > 8 kcal/mol but pose passes PoseBusters geometric.

**Fix:** Re-dock with `exhaustiveness=32`; or post-relax pose with UFF inside pocket.

## Reconciliation: PoseBusters vs RMSD

| RMSD <= 2A | PB-valid | Action |
|------------|----------|--------|
| Yes | Yes | High confidence; use for SAR |
| Yes | No | RMSD-good but unphysical; re-relax |
| No | Yes | Physical but wrong binding mode; consider ensemble |
| No | No | Both wrong; reject |

In DiffDock benchmarks, ~25% of poses are RMSD <= 2A AND PB-valid; ~25% RMSD <= 2A but NOT PB-valid (the "false success" rate of AI methods).

## Integration into VS Pipeline

```python
import pandas as pd
from posebusters import PoseBusters

def pose_qc_pipeline(docked_sdfs, receptor_pdb):
    bust = PoseBusters(config='dock')
    all_results = []
    for sdf in docked_sdfs:
        r = bust.bust(mol_pred=sdf, mol_cond=receptor_pdb)
        r['pb_valid'] = r.iloc[:, 4:].all(axis=1)
        r['source'] = sdf
        all_results.append(r)
    df = pd.concat(all_results)

    df['rank'] = df.groupby('source')['pb_valid'].cumsum()
    valid_top = df[df['pb_valid']].groupby('source').head(1)
    return valid_top
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| PoseBusters reports None | Sanitize failed on input pose | Standardize before docking; check SDF format |
| RMSD not computed | No reference provided | Pass `mol_true` parameter |
| All checks pass for invalid pose | Wrong receptor file format | Use PDB with hydrogens; PDBQT may not work |
| vdW overlap false positive on covalent | Covalent bond counted as clash | Use covalent docking-specific validation |
| Strain calculation slow | Too many reference conformers | Reduce `n_ref` to 5-10 |
| Posebusters config error | Wrong config name | Options: `redock`, `dock`, `mol` |
| posecheck unavailable | Different tool, similar purpose | `pip install posecheck` for alternative |

## References

- Buttenschoen et al., *Chem. Sci.* 15:3130 (2024) -- PoseBusters benchmark and tool.
- Boström, *J. Comput.-Aided Mol. Des.* 12:383 (1998) -- bioactive ligand strain.
- Hawkins, *J. Chem. Inf. Model.* 57:1747 (2017) -- conformer ensemble best practices.
- Krasoulis et al., *J. Cheminformatics* 14:24 (2022) -- DENVIS pose quality assessment.
- Cole et al., *Acta Crystallogr. D* 64:144 (2008) -- ligand strain in PDB structures.

## Related Skills

- chemoinformatics/virtual-screening - Source of poses to validate
- chemoinformatics/ml-docking-rescoring - DiffDock, EquiBind, TANKBind validation
- chemoinformatics/molecular-io - SDF format handling
- chemoinformatics/conformer-generation - Generate reference conformer ensemble for strain
- chemoinformatics/free-energy-calculations - PoseBusters-valid poses for FEP input
- chemoinformatics/covalent-design - Covalent pose validation
