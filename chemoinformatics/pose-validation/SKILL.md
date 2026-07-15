---
name: bio-pose-validation
description: Validates docked / generated protein-ligand poses using PoseBusters physical-validity tests, strain energy quantification, geometric checks (planarity, vdW overlap, bond/angle distortion), and pose-energy reasonableness. Use when QC-ing docking results, comparing classical vs ML docking outputs, or filtering pose lists before SAR analysis.
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

Test docked or AI-generated protein-ligand poses for physical plausibility. PoseBusters (Buttenschoen et al. 2024) provides geometric, chemical, and energetic checks that flag implausible poses, including non-planar aromatic rings, van der Waals clashes, broken bonds, altered stereochemistry, and unfavorable internal energies. On the Astex Diverse Set, DiffDock achieved 72% RMSD success but only 47% combined RMSD-and-PB-valid success; the size of this gap is dataset- and method-dependent. PB-valid status complements RMSD for downstream SAR, FEP setup, or generative-model training.

For docking, see `chemoinformatics/virtual-screening`. For ML docking specifically, see `chemoinformatics/ml-docking-rescoring`.

## PoseBusters Test Suite

PoseBusters runs ~20 individual checks grouped into:

The thresholds below are the benchmark criteria reported by Buttenschoen et al. (2024). Installed PoseBusters defaults may differ by version and configuration, so record the package version and resolved configuration.

| Check group | What it tests | 2024 benchmark criterion |
|-------------|---------------|-----------|
| Sanity | Ligand chemical sanity | RDKit sanitization passes |
| Bond lengths | Bond lengths within reference | 0.75–1.25 times RDKit distance-geometry bounds |
| Bond angles | 1–3 distances within reference | 0.75–1.25 times RDKit distance-geometry bounds |
| Internal steric | No intra-ligand clash | Pair distance > 0.70 times the RDKit lower bound |
| Aromatic ring planarity | Aromatic rings planar | Maximum deviation from fitted plane <= 0.25 Å |
| Double-bond stereo | Z/E preserved | Match input SMILES |
| Internal energy | Energy relative to generated conformers | UFF energy ratio <= 100 versus the mean of 50 generated, relaxed conformers |
| Volume overlap | vdW overlap with protein | < 7.5% of ligand vdW volume |
| Minimum distance | No severe protein-ligand clash | Distance >= 0.75 times the sum of vdW radii |
| Chirality | R/S preserved from input | Match input SMILES |

A pose passing ALL tests is "PB-valid". Combined PB-valid + RMSD <= 2 Å is the modern criterion.

## When to Apply PoseBusters

| Workflow | PoseBusters use | Action |
|----------|-----------------|--------|
| Self-docking (validating method) | Required | Compare PB-valid + RMSD <= 2A |
| Cross-docking | Required | PB-valid + RMSD <= 2A; account for protein flexibility |
| Virtual screening top hits | Required | Filter to PB-valid before MM/GBSA / FEP |
| AI docking (DiffDock, etc.) | Required for a fair benchmark | Report the dataset-specific PB-valid and combined success rates |
| Generated ligand poses | Recommended | Measure chemical and geometric validity rather than assuming it |
| Boltz-2 / AlphaFold3 ligand poses | Recommended | Benchmark validity on the relevant complexes; do not infer a failure frequency from DiffDock |
| Production FEP setup | Required | Inspect pose validity and ligand strain before system preparation |

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

Common configurations and their included checks are:

| Config | Includes | When to use |
|--------|----------|-------------|
| `redock` | All checks + RMSD vs reference + protein vdW overlap | Self-docking benchmarks, retrospective validation |
| `dock` | All checks except RMSD reference | Blind docking, prospective virtual screening |
| `mol` | Intra-ligand only (sanity, bonds, angles, rings, stereo, energy) | Conformer QC; no protein context |

PoseBusters also ships additional and faster configurations in some releases. Treat the table as a workflow guide, not an exhaustive registry, and inspect the configurations available in the installed version.

Output: a DataFrame with one row per pose, metadata columns, and boolean pass/fail columns for the checks enabled by the selected configuration. Reference-dependent fields such as RMSD and the exact check-column names vary by configuration and version; inspect `results.columns` rather than relying on a fixed exhaustive list.

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

check_cols = [
    col for col in results.select_dtypes(include='bool').columns
    if not col.lower().startswith('rmsd')
]
results['pb_valid'] = results[check_cols].all(axis=1)
valid = results[results['pb_valid']]
print(f'{len(valid)} / {len(results)} poses are PB-valid')
```

## Strain Energy Quantification

Beyond binary PB-valid, quantitative strain energy distinguishes "marginal" from "egregious" poses.

**Goal:** Quantify how far each docked pose is from its lowest-energy free conformer in MMFF94 energy units.

**Approach:** Generate a reference conformer ensemble (ETKDGv3 + MMFF94), make the docked and reference molecules chemically consistent by adding explicit hydrogens to both, relax only the added docked-pose hydrogens while fixing all heavy atoms, take the lowest sampled reference energy as baseline, and report `docked_energy - min_ref_energy` as a relative strain diagnostic. This is not a rigorous solution-phase conformational free energy.

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
        if ref is None:
            strains.append({'strain': None, 'note': 'reference_parse_failed'})
            continue
        ref = Chem.AddHs(ref)
        props_ref = AllChem.MMFFGetMoleculeProperties(ref)
        if props_ref is None:
            strains.append({'strain': None, 'note': 'no_reference_mmff_parameters'})
            continue
        conf_ids = list(AllChem.EmbedMultipleConfs(
            ref, numConfs=n_ref, params=AllChem.ETKDGv3()
        ))
        if not conf_ids:
            strains.append({'strain': None, 'note': 'reference_embedding_failed'})
            continue
        AllChem.MMFFOptimizeMoleculeConfs(ref)

        ref_energies = []
        for c in conf_ids:
            ff = AllChem.MMFFGetMoleculeForceField(
                ref, props_ref, confId=c
            )
            if ff is not None:
                ref_energies.append(ff.CalcEnergy())
        if not ref_energies:
            strains.append({'strain': None, 'note': 'reference_force_field_failed'})
            continue
        min_ref = min(ref_energies)

        # MMFF energies are comparable only for the same explicit atom system.
        # Add any missing H coordinates, then relax H atoms while preserving the
        # docked heavy-atom pose.
        docked_h = Chem.AddHs(Chem.Mol(docked), addCoords=True)
        if docked_h.GetNumAtoms() != ref.GetNumAtoms():
            strains.append({'strain': None, 'note': 'atom_system_mismatch'})
            continue
        props_docked = AllChem.MMFFGetMoleculeProperties(docked_h)
        docked_ff = AllChem.MMFFGetMoleculeForceField(
            docked_h, props_docked
        ) if props_docked is not None else None
        if docked_ff is not None:
            for atom in docked_h.GetAtoms():
                if atom.GetAtomicNum() != 1:
                    docked_ff.AddFixedPoint(atom.GetIdx())
            docked_ff.Minimize(maxIts=200)
        docked_e = docked_ff.CalcEnergy() if docked_ff else None

        strains.append({
            'min_ref_energy': min_ref,
            'docked_energy': docked_e,
            'strain': docked_e - min_ref if docked_e is not None else None,
            'note': 'ok' if docked_e is not None else 'docked_force_field_failed',
        })
    return strains
```

Interpret relative MMFF strain in the context of ligand chemistry, conformer-sampling coverage, and force-field support. Boström et al. (1998) found a conformational energy penalty of no more than 3 kcal/mol for about 70% of 33 protein-bound ligands; that result does not establish a universal acceptance cutoff. Treat unusually high values as a prompt for inspection or use a project-defined threshold validated for the series.

## vdW Overlap with Protein

The 2024 benchmark criterion limits protein-ligand overlap to 7.5% of the ligand vdW volume, using protein radii scaled by 0.8. PoseBusters' `bust(...)` computes this check; do not substitute an unvalidated pairwise-distance sketch for its volume calculation.

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

## Model-Specific Failure Diagnosis

Do not assign a mechanism from the model name or a failed PoseBusters column alone. For DiffDock-L, EquiBind, TANKBind, Boltz, AlphaFold3, or another pose generator, report the observed failed checks on the evaluated dataset, inspect the structures, and compare against the method's documented constraints. A chirality, planarity, bond-geometry, or clash failure may justify filtering or a validated constrained-relaxation protocol, but relaxation must be checked for displacement of the binding mode.

### High strain after Vina docking

**Trigger:** Highly constrained pocket; flexible ligand.

**Symptom:** Relative strain is an outlier for the chemical series even though the pose passes the enabled geometric checks.

**Fix:** Inspect conformer-sampling coverage and force-field support. Compare additional docking or constrained-relaxation settings under a project-validated protocol rather than applying a universal strain or exhaustiveness cutoff.

## Reconciliation: PoseBusters vs RMSD

| RMSD <= 2A | PB-valid | Action |
|------------|----------|--------|
| Yes | Yes | Physically plausible and close to the reference; still validate suitability for the downstream task |
| Yes | No | Close to the reference but fails an enabled plausibility check; inspect the failure and any validated relaxation |
| No | Yes | Physically plausible but different from the reference; investigate alignment, protein state, and alternative binding modes |
| No | No | Different from the reference and fails an enabled plausibility check; inspect both causes before deciding whether to reject |

On the Astex Diverse Set reported by Buttenschoen et al. (2024), DiffDock's top-pose success fell from 72% by RMSD <= 2 Å alone to 47% when PB-validity was also required: a 25-percentage-point gap. Do not generalize that result to a fixed failure rate on other datasets.

## Integration into VS Pipeline

```python
import pandas as pd
from posebusters import PoseBusters

def pose_qc_pipeline(docked_sdfs, receptor_pdb):
    bust = PoseBusters(config='dock')
    all_results = []
    for sdf in docked_sdfs:
        r = bust.bust(mol_pred=sdf, mol_cond=receptor_pdb)
        check_cols = [
            col for col in r.select_dtypes(include='bool').columns
            if not col.lower().startswith('rmsd')
        ]
        r['pb_valid'] = r[check_cols].all(axis=1)
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
| Rows or expected checks are missing | Input loading failed or the selected configuration omits those checks | Inspect the returned DataFrame, loading-status columns, input format, and installed configuration |
| RMSD not computed | No reference provided | Pass `mol_true` parameter |
| All checks pass for invalid pose | Wrong receptor file format | Use PDB with hydrogens; PDBQT may not work |
| vdW overlap false positive on covalent | Covalent bond counted as clash | Use covalent docking-specific validation |
| Strain calculation slow | Too many reference conformers | Reduce `n_ref` to 5-10 |
| PoseBusters config error | Wrong or version-incompatible config name | Inspect the installed configuration registry; `redock`, `dock`, and `mol` are common configurations |
| posecheck unavailable | Different tool, similar purpose | `pip install posecheck` for alternative |

## References

- Buttenschoen M, Morris GM, Deane CM. "PoseBusters: AI-based docking methods fail to generate physically valid poses or generalise to novel sequences." *Chem. Sci.* 15:3130–3139 (2024). DOI: 10.1039/D3SC04185A.
- Boström J, Norrby PO, Liljefors T. "Conformational energy penalties of protein-bound ligands." *J. Comput.-Aided Mol. Des.* 12:383–396 (1998). DOI: 10.1023/A:1008007507641.
- Corso G et al. "DiffDock: Diffusion Steps, Twists, and Turns for Molecular Docking." *ICLR* (2023). OpenReview: https://openreview.net/forum?id=kKF8_K-mBbS.
- Stärk H et al. "EquiBind: Geometric Deep Learning for Drug Binding Structure Prediction." *PMLR* 162:20503–20521 (2022). https://proceedings.mlr.press/v162/stark22b.html.
- Lu W et al. "TankBind: Trigonometry-Aware Neural NetworKs for Drug-Protein Binding Structure Prediction." *NeurIPS* 35 (2022). Official repository: https://github.com/luwei0917/TankBind.
- Abramson J et al. "Accurate structure prediction of biomolecular interactions with AlphaFold 3." *Nature* 630:493–500 (2024). DOI: 10.1038/s41586-024-07487-w.
- Boltz official repository and documentation: https://github.com/jwohlwend/boltz.
- PoseBusters documentation, Python API: https://posebusters.readthedocs.io/en/latest/api.html.

## Related Skills

- chemoinformatics/virtual-screening - Source of poses to validate
- chemoinformatics/ml-docking-rescoring - DiffDock, EquiBind, TANKBind validation
- chemoinformatics/molecular-io - SDF format handling
- chemoinformatics/conformer-generation - Generate reference conformer ensemble for strain
- chemoinformatics/free-energy-calculations - PoseBusters-valid poses for FEP input
- chemoinformatics/covalent-design - Covalent pose validation
