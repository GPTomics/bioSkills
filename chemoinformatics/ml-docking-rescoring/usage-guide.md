# ML Docking and Rescoring Usage Guide

## Overview

Modern ML-based protein-ligand pose prediction and scoring. DiffDock-L (diffusion), Boltz-1 / Boltz-2 (foundation model with affinity), Chai-1, AlphaFold3 ligand, EquiBind, TANKBind, NeuralPLexer. Modern best practice: hybrid pipeline (DiffDock-L pose + GNINA rescore + PoseBusters QC).

## Prerequisites

```bash
pip install rdkit posebusters gnina
# Boltz weights: from official repo
# Chai-1: pip install chai-lab
# DiffDock-L: GitHub installation
```

## Quick Start

Tell the AI agent what to do:
- "Dock my ligand with DiffDock-L; validate with PoseBusters"
- "Predict affinity for 1000 candidates with Boltz-2 affinity module"
- "Run AlphaFold3 ligand prediction for novel target"
- "Hybrid: DiffDock pose + GNINA rescore + PoseBusters filter"

## Example Prompts

### Hybrid VS pipeline
> "For 100 ligands in lib.csv against receptor.pdb: DiffDock-L sample 5 poses each -> GNINA CNN rescore -> PoseBusters filter. Output top 10 PB-valid + RMSD-validated."

### Boltz-2 affinity triage
> "Boltz-2 affinity prediction for 1000 SMILES against PDB 4XYZ. Rank by predicted affinity. Output top 50 for FEP follow-up."

### Cross-docking
> "Cross-dock 50 ligands against AlphaFold-predicted holo. Use DiffDock-L (best for cross-dock); validate physical plausibility."

### Boltz-1 multimer
> "Predict ternary complex (target + PROTAC + E3) with Boltz-1. Use chain-chain restraints."

## What the Agent Will Do

1. Set up DiffDock-L / Boltz / Chai-1 / AlphaFold3 input (PDB + SMILES).
2. Generate 5-40 poses per ligand (model-dependent).
3. Rescore poses with GNINA CNN scoring.
4. Run PoseBusters per pose; filter PB-valid.
5. Rank by combined score (DiffDock confidence + GNINA + Boltz-2 affinity).
6. Output PB-valid + ranked candidates.

## Tips

- DiffDock-L produces ~50% PB-invalid poses; PoseBusters validation is mandatory.
- Boltz-2 affinity: Pearson 0.66 on FEP benchmark; 1000x faster.
- Chai-1 / AlphaFold3 best for novel target with limited co-crystal data.
- For ultralarge libraries, use classical Vina pre-filter + ML rescore on top.
- Always validate with experimental binding when possible.

## Related Skills

- chemoinformatics/virtual-screening - Classical docking foundation
- chemoinformatics/pose-validation - PoseBusters QC (mandatory after ML docking)
- chemoinformatics/free-energy-calculations - Boltz-2 alternative
- chemoinformatics/molecular-io - Format conversion
- chemoinformatics/conformer-generation - Pre-conformer for some ML tools
- structural-biology/modern-structure-prediction - Protein structure prediction
- structural-biology/structure-io - PDB / mmCIF
