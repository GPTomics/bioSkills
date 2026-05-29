# Free Energy Calculations Usage Guide

## Overview

Perform alchemical free-energy calculations: RBFE (relative binding) and ABFE (absolute binding) using OpenFE (open-source) or FEP+ (commercial). Compute thermodynamic-cycle-validated affinity differences with 1-2 kcal/mol RMSE.

## Prerequisites

```bash
conda install -c conda-forge openff-toolkit openmm gromacs
pip install openfe alchemlyb pymbar
```

## Quick Start

Tell the AI agent what to do:
- "Run OpenFE RBFE between lig1 and lig2 in this receptor"
- "Compute ABFE for a single ligand bound to my receptor"
- "Run MM/GBSA endpoint scoring on docked poses"
- "Check FEP cycle closure for my 5-ligand perturbation graph"
- "Boltz-2 affinity first-pass on 100 candidates"

## Example Prompts

### Lead optimization RBFE
> "Run OpenFE RBFE on a 6-ligand series for kinase X. Use 12 lambda windows, 5 ns/window. Report delta-delta-G with cycle closure error. Output Forcefield: SAGE 2.1.0."

### Quick MM/GBSA screen
> "Run MM/GBSA on the top 100 Vina-docked poses for receptor.pdb. Rank by delta-G_binding. Output to results.csv with sd."

### Boltz-2 affinity triage
> "Boltz-2 affinity prediction on 1000 SMILES against receptor.pdb. Rank by predicted affinity. Output top 20 for FEP follow-up."

### ABFE for single ligand
> "Compute OpenFE ABFE for ligand.sdf bound to receptor.pdb. Use Boresch-style restraints on rigid scaffold atoms; 5 ns charging, 12 ns vdW, 7 ns restraint windows."

## What the Agent Will Do

1. Set up perturbation: ligand1 -> ligand2 (RBFE) or ligand -> uncoupled (ABFE).
2. Build hybrid topology with atom mapping (LOMAP for RBFE).
3. Generate lambda window schedule (12-20 windows typical).
4. Equilibrate each window, then production MD (5-20 ns/window).
5. Analyze with MBAR/BAR via alchemlyb.
6. Compute cycle closure error; report delta-delta-G ± SD.
7. Validate convergence: replicates should agree within 0.5 kcal/mol.

## Tips

- OpenFE is open-source equivalent to FEP+; SAGE 2.1.0 default forcefield.
- Always run REST2 enhanced sampling for buried pockets.
- Cycle closure > 1 kcal/mol RMS indicates insufficient sampling.
- ABFE is ~3x cost of RBFE; use selectively.
- MM/GBSA is fast first-pass (3-5 kcal/mol RMSE), not lead-optimization standard.
- Boltz-2 affinity approaches FEP accuracy (Pearson 0.66) at 1000x speed.

## Related Skills

- chemoinformatics/virtual-screening - Source docking poses for FEP input
- chemoinformatics/pose-validation - PoseBusters before FEP setup
- chemoinformatics/conformer-generation - Ligand 3D for FEP
- chemoinformatics/molecular-standardization - Standardize before FEP
- chemoinformatics/ml-docking-rescoring - Boltz-2 affinity alternative
- chemoinformatics/qsar-modeling - Surrogate models for screening
