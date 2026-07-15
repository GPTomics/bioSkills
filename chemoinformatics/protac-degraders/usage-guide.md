# PROTAC and Bivalent Degrader Design Usage Guide

## Overview

Design PROTACs (bivalent molecules recruiting E3 ligase to target for proteasomal degradation). Balance target-ligand, E3-ligand, linker geometry, cooperativity, and cell permeability. Cover ternary complex prediction (PRosettaC, DeepTernary, AlphaFold3), cooperativity, hook effect, and DC50/Dmax.

## Prerequisites

```bash
pip install rdkit
# PRosettaC: web service (prosettac.weizmann.ac.il)
# DeepTernary: GitHub installation
# AlphaFold3: use AlphaFold Server or a licensed local installation; arbitrary distance restraints are not supported
```

## Quick Start

Tell the AI agent what to do:
- "Design PROTAC for kinase X target using CRBN E3"
- "Enumerate linkers between target-ligand and VHL ligand"
- "Predict ternary complex for kinase + PROTAC + CRBN"
- "Optimize linker length to maximize cooperativity"

## Example Prompts

### CRBN PROTAC design
> "For target kinase X (PDB 5XYZ), design an exploratory PROTAC series using a justified E3 recruiter. Vary linker composition and geometry around candidates compatible with the binary structures. Output connected, synthesis-review-ready SMILES with structural scores reported as hypotheses."

### Ternary complex prediction
> "Predict the ternary structure for this target-ligand-E3 PROTAC complex using PRosettaC. Report structural and interface scores; state that cooperativity alpha requires a binding experiment."

### Linker optimization
> "Given target-ligand and E3-ligand exit vectors in a shared structural hypothesis, suggest a small linker series spanning rigid and flexible chemistries. Explain how conformer feasibility will be checked rather than inferring atom count from distance alone."

### VHL alternative design
> "Switch from CRBN to VHL E3. Adjust linker to maintain ternary geometry. Re-predict ternary complex."

## What the Agent Will Do

1. Define target ligand (from co-crystal or docked) and E3 ligand (pomalidomide / VHL ligand).
2. Compute distance between attachment points on each ligand.
3. Enumerate a linker series around geometries supported by the binary structures.
4. Combine target-linker-E3 SMILES.
5. Predict ternary complex via PRosettaC (or DeepTernary).
6. Score by linker geometry and structural/interface metrics; report experimental cooperativity separately when available.

## Tips

- CRBN and VHL have extensively published recruiter series, but choose the E3 using target/E3 geometry, expression in the intended system, ligand availability, and known liabilities.
- Optimal linker length is target-specific; use structural hypotheses to choose an exploratory series rather than a universal range.
- Positive cooperativity may be beneficial, but alpha is assay- and system-dependent and must be measured experimentally.
- Hook effect at high concentration; bell-shaped dose-response.
- Permeability can be a bottleneck; compare measured permeability and cellular activity across the series rather than applying a universal MW or TPSA cutoff.
- The linker-enumeration example requires each explicit dummy attachment to use a single bond; reject other dummy-bond orders unless the intended connection chemistry has a separately validated bond-order rule.

## Related Skills

- chemoinformatics/molecular-io - Parse ligand SMILES
- chemoinformatics/reaction-enumeration - Linker enumeration
- chemoinformatics/generative-design - REINVENT linker mode
- chemoinformatics/conformer-generation - Ternary conformer sampling
- chemoinformatics/virtual-screening - Validate target ligand binding
- chemoinformatics/free-energy-calculations - Ternary ABFE
- chemoinformatics/admet-prediction - PROTAC ADMET specifics
- structural-biology/structure-io - PDB / mmCIF for ternary complex
