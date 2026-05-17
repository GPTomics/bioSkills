# Pose Validation Usage Guide

## Overview

Validate docked or AI-generated protein-ligand poses for physical plausibility using PoseBusters (the modern standard). Quantitatively measure ligand strain, geometric distortion, vdW overlap, and stereochemistry preservation. Essential for filtering AI docking outputs (DiffDock, EquiBind), where ~50% of poses fail physical-validity tests.

## Prerequisites

```bash
pip install posebusters rdkit pandas
```

## Quick Start

Tell the AI agent what to do:
- "Run PoseBusters on my docked SDF and filter to PB-valid poses"
- "Compute ligand strain energy for each docked pose"
- "Compare DiffDock vs GNINA pose validity on the PoseBusters benchmark"
- "Validate poses before FEP setup: strain < 5 kcal/mol, PB-valid, RMSD < 2 A"
- "Identify chirality-inverted poses in my DiffDock output"

## Example Prompts

### Standard pose QC
> "Run PoseBusters dock config on docked_poses.sdf with receptor.pdb. Output a table per pose with each check pass/fail. Filter to PB-valid; rank by Vina score."

### Strain energy quantification
> "For each docked pose, compute the MMFF94 strain energy vs the minimized free conformer. Flag poses with strain > 5 kcal/mol as potential errors."

### AI docking validation
> "Run DiffDock-L on 100 ligands. For each, run PoseBusters; report PB-valid rate. Compare to GNINA classical docking on same compounds."

### FEP input prep
> "Filter docked poses to those passing PoseBusters AND strain < 5 kcal/mol. Output as FEP-ready SDF."

## What the Agent Will Do

1. Read input SDF of docked poses and PDB of receptor.
2. Run PoseBusters bust() for each pose against receptor.
3. Combine all PoseBusters checks into PB-valid boolean.
4. Optionally compute strain energy per pose.
5. Filter to PB-valid; report which checks failed for invalid poses.
6. Rank surviving poses by docking score or affinity prediction.

## Tips

- Always validate AI docking output (DiffDock-L, EquiBind, TANKBind) -- 50% fail rate typical.
- PB-valid + RMSD <= 2 Å is the modern gold standard.
- Strain energy: < 5 kcal/mol acceptable; > 10 kcal/mol implausible.
- For covalent docking, exclude vdW overlap check.
- Use PoseBusters as filter, not absolute reject; sometimes valid poses fail on edge cases.

## Related Skills

- chemoinformatics/virtual-screening - Source poses from classical docking
- chemoinformatics/ml-docking-rescoring - DiffDock + GNINA + PoseBusters hybrid
- chemoinformatics/molecular-io - SDF parsing
- chemoinformatics/conformer-generation - Reference conformers for strain
- chemoinformatics/free-energy-calculations - PB-valid input for FEP
- chemoinformatics/covalent-design - Covalent pose validation
