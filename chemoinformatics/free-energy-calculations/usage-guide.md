# Free Energy Calculations Usage Guide

## Overview

Perform alchemical free-energy calculations using OpenFE or a commercial workflow. Accuracy is system- and protocol-dependent, so validate mappings, sampling, network consistency, and prospective performance rather than promising a universal RMSE.

## Prerequisites

```bash
mamba create -n openfe -c conda-forge openfe=1.7 alchemlyb pymbar
conda activate openfe
```

Install OpenFE from conda-forge as documented by the project. The unrelated package named `openfe` on PyPI is not the Open Free Energy toolkit.

## Quick Start

Tell the AI agent what to do:
- "Run OpenFE RBFE between lig1 and lig2 in this receptor"
- "Compute ABFE for a single ligand bound to my receptor"
- "Run MM/GBSA endpoint scoring on docked poses"
- "Check FEP cycle closure for my 5-ligand perturbation graph"
- "Boltz-2 affinity first-pass on 100 candidates"

## Example Prompts

### Lead optimization RBFE
> "Run OpenFE RBFE on a congeneric series for kinase X. Start from versioned protocol defaults, inspect mappings, and extend sampling from overlap, replicate, and closure-residual diagnostics."

### Quick MM/GBSA screen
> "Run MM/GBSA on the top 100 Vina-docked poses for receptor.pdb. Rank by delta-G_binding. Output to results.csv with sd."

### Boltz-2 affinity triage
> "Boltz-2 affinity prediction on 1000 SMILES against receptor.pdb. Rank by predicted affinity. Output top 20 for FEP follow-up."

### ABFE for single ligand
> "Compute OpenFE ABFE for ligand.sdf bound to receptor.pdb using the documented restraint workflow. Record all settings and extend each leg from convergence diagnostics."

## What the Agent Will Do

1. Set up perturbation: ligand1 -> ligand2 (RBFE) or ligand -> uncoupled (ABFE).
2. Generate and inspect atom mappings (Kartograf is the OpenFE 1.7 CLI default; LOMAP is supported).
3. Start from the release-matched lambda schedule and modify it only from overlap and protocol diagnostics.
4. Equilibrate and run production sampling, extending from replicate and time-series convergence evidence.
5. Analyze with MBAR/BAR via alchemlyb.
6. Compute cycle closure error; report delta-delta-G ± SD.
7. Validate convergence using replicate agreement, overlap, exchange, and signed closure residuals with propagated uncertainty.

## Tips

- OpenFE 1.7's versioned documentation shows OpenFF 2.1.1 for ligands; inspect and record the actual installed protocol settings.
- OpenFE's default RBFE Hamiltonian replica exchange is not REST2. Use enhanced sampling only when the selected engine and protocol document it.
- Interpret signed cycle-closure residuals with edge uncertainties, shared-edge correlations, and replicate diagnostics; RMS requires multiple specified cycles.
- Estimate ABFE/RBFE cost from the explicit protocol rather than a fixed multiplier.
- Use MM/GBSA only where a matched benchmark supports the intended ranking decision.
- The Boltz-2 preprint reports Pearson 0.66 on its FEP benchmark subset and at least 1000-fold lower computational cost; treat this as benchmark-specific and confirm top candidates with FEP.

## Related Skills

- chemoinformatics/virtual-screening - Source docking poses for FEP input
- chemoinformatics/pose-validation - PoseBusters before FEP setup
- chemoinformatics/conformer-generation - Ligand 3D for FEP
- chemoinformatics/molecular-standardization - Standardize before FEP
- chemoinformatics/ml-docking-rescoring - Boltz-2 affinity alternative
- chemoinformatics/qsar-modeling - Surrogate models for screening
