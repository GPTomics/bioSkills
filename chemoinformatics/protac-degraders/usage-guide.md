# PROTAC and Bivalent Degrader Design Usage Guide

## Overview

Design PROTACs (bivalent molecules recruiting E3 ligase to target for proteasomal degradation). Balance target-ligand, E3-ligand, linker geometry, cooperativity, and cell permeability. Cover ternary complex prediction (PRosettaC, DeepTernary, AlphaFold3), cooperativity, hook effect, and DC50/Dmax.

## Prerequisites

```bash
pip install rdkit
# PRosettaC: web service (prosettac.weizmann.ac.il)
# DeepTernary: GitHub installation
# AlphaFold3: API access required
```

## Quick Start

Tell the AI agent what to do:
- "Design PROTAC for kinase X target using CRBN E3"
- "Enumerate linkers between target-ligand and VHL ligand"
- "Predict ternary complex for kinase + PROTAC + CRBN"
- "Optimize linker length to maximize cooperativity"

## Example Prompts

### CRBN PROTAC design
> "For target kinase X (PDB 5XYZ), design 30 PROTACs using CRBN E3 (pomalidomide). Vary linker length 10-20 atoms with PEG / piperazine / triazole. Output SMILES sorted by predicted ternary stability."

### Ternary complex prediction
> "Predict ternary complex for given target-ligand-E3 PROTAC complex. Use PRosettaC. Report cooperativity (alpha) estimate."

### Linker optimization
> "Given target ligand exit vector and CRBN ligand entry vector, suggest 5 linkers of optimal length. Use rigid (piperazine) and flexible (PEG) variants."

### VHL alternative design
> "Switch from CRBN to VHL E3. Adjust linker to maintain ternary geometry. Re-predict ternary complex."

## What the Agent Will Do

1. Define target ligand (from co-crystal or docked) and E3 ligand (pomalidomide / VHL ligand).
2. Compute distance between attachment points on each ligand.
3. Enumerate linkers within target distance range (e.g., 12-18 atoms).
4. Combine target-linker-E3 SMILES.
5. Predict ternary complex via PRosettaC (or DeepTernary).
6. Score by linker geometry, ternary stability, Lipinski compliance.

## Tips

- CRBN is the most-developed E3 (broad tissue, thalidomide-derived ligands).
- VHL is the second-most-developed (more selective, tissue-restricted).
- Optimal linker length is target-specific; typically 12-20 atoms.
- Positive cooperativity (alpha > 2) is gold standard.
- Hook effect at high concentration; bell-shaped dose-response.
- Permeability is bottleneck; PROTAC 800-1200 Da, TPSA < 130 typical.

## Related Skills

- chemoinformatics/molecular-io - Parse ligand SMILES
- chemoinformatics/reaction-enumeration - Linker enumeration
- chemoinformatics/generative-design - REINVENT linker mode
- chemoinformatics/conformer-generation - Ternary conformer sampling
- chemoinformatics/virtual-screening - Validate target ligand binding
- chemoinformatics/free-energy-calculations - Ternary ABFE
- chemoinformatics/admet-prediction - PROTAC ADMET specifics
- structural-biology/structure-io - PDB / mmCIF for ternary complex
