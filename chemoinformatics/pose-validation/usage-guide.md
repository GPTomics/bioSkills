# Pose Validation Usage Guide

## Overview

Validate docked or AI-generated protein-ligand poses for physical plausibility using PoseBusters. Quantitatively measure ligand strain, geometric distortion, vdW overlap, and stereochemistry preservation. Report dataset-specific PB-valid rates rather than assuming a fixed failure rate for an AI docking method.

## Prerequisites

```bash
pip install posebusters rdkit pandas
```

## Quick Start

Tell the AI agent what to do:
- "Run PoseBusters on my docked SDF and filter to PB-valid poses"
- "Compute ligand strain energy for each docked pose"
- "Compare DiffDock vs GNINA pose validity on the PoseBusters benchmark"
- "Validate poses before FEP setup: PB-valid, RMSD < 2 Å when a reference is available, and report relative MMFF strain"
- "Identify chirality-inverted poses in my DiffDock output"

## Example Prompts

### Standard pose QC
> "Run PoseBusters dock config on docked_poses.sdf with receptor.pdb. Output a table per pose with each check pass/fail. Filter to PB-valid; rank by Vina score."

### Strain energy quantification
> "For each docked pose, compute relative MMFF94 strain versus the lowest sampled reference conformer. Report high-strain outliers for inspection without imposing an unvalidated universal cutoff."

### AI docking validation
> "Run DiffDock-L on 100 ligands. For each, run PoseBusters; report PB-valid rate. Compare to GNINA classical docking on same compounds."

### FEP input prep
> "Filter docked poses to those passing PoseBusters. Report relative MMFF strain and apply only the project-defined cutoff validated for this chemical series before FEP setup."

## What the Agent Will Do

1. Read input SDF of docked poses and PDB of receptor.
2. Run PoseBusters bust() for each pose against receptor.
3. Combine all PoseBusters checks into PB-valid boolean.
4. Optionally compute strain energy per pose.
5. Filter to PB-valid; report which checks failed for invalid poses.
6. Rank surviving poses by docking score or affinity prediction.

## Tips

- Always validate AI docking output and report the observed failure rate for the dataset.
- In redocking benchmarks, report PB-validity together with RMSD rather than RMSD alone.
- Relative MMFF strain is a diagnostic whose interpretation depends on force-field support and conformer sampling. Compare identical explicit atom systems: add hydrogens consistently to the docked and reference molecules, relaxing only docked-pose hydrogens while fixing heavy atoms. Use a validated project-specific cutoff if filtering on it.
- For covalent docking, exclude vdW overlap check.
- Use PoseBusters as filter, not absolute reject; sometimes valid poses fail on edge cases.

## Related Skills

- chemoinformatics/virtual-screening - Source poses from classical docking
- chemoinformatics/ml-docking-rescoring - DiffDock + GNINA + PoseBusters hybrid
- chemoinformatics/molecular-io - SDF parsing
- chemoinformatics/conformer-generation - Reference conformers for strain
- chemoinformatics/free-energy-calculations - PB-valid input for FEP
- chemoinformatics/covalent-design - Covalent pose validation
