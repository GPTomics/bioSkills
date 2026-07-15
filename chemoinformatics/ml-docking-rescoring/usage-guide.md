# ML Docking and Rescoring Usage Guide

## Overview

Modern ML-based protein-ligand pose prediction and scoring. DiffDock-L (diffusion), Boltz-1 / Boltz-2 (foundation model with affinity), Chai-1, AlphaFold3 ligand, EquiBind, TANKBind, and NeuralPLexer are candidate samplers or predictors. An auditable hybrid evaluation can keep DiffDock-L poses, GNINA rescoring, and PoseBusters QC as separate evidence streams.

## Prerequisites

```bash
pip install rdkit posebusters
# Install the GNINA executable from its documented binary/container distribution.
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
> "For 100 ligands in lib.csv against receptor.pdb: DiffDock-L sample 5 poses each -> GNINA CNN rescore -> PoseBusters filter. Output the top 10 PB-valid poses; calculate RMSD only if I provide a reference pose."

### Boltz-2 affinity triage
> "Boltz-2 affinity prediction for 1000 SMILES against PDB 4XYZ. Rank by predicted affinity. Output top 50 for FEP follow-up."

### Cross-docking
> "Cross-dock 50 ligands against an AlphaFold-predicted receptor with DiffDock-L, then rescore and validate physical plausibility."

### Boltz-1 multimer
> "Predict ternary complex (target + PROTAC + E3) with Boltz-1. Use chain-chain restraints."

## What the Agent Will Do

1. Set up DiffDock-L, Boltz, Chai-1, or AlphaFold 3 input using that release's official schema.
2. Generate 5-40 poses per ligand (model-dependent).
3. Rescore poses with GNINA CNN scoring.
4. Run PoseBusters per pose; filter PB-valid.
5. Keep model confidence, GNINA score, optional Boltz-2 affinity, and physical-validity results separate; use consensus ranks and inspect disagreements.
6. Output PB-valid + ranked candidates.

## Tips

- Pose confidence does not establish physical validity; apply PoseBusters checks to generated poses.
- Boltz-2 reports affinity performance approaching FEP on evaluated benchmarks at at least 1,000-fold lower computational cost; do not transfer a single benchmark correlation to new targets. Use `affinity_probability_binary` for hit discovery and `affinity_pred_value` for hit-to-lead or lead optimization, with target-relevant validation.
- Chai-1 and AlphaFold 3 are options when a complex structure must be predicted; compare confidence and validate the ligand pose independently.
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
