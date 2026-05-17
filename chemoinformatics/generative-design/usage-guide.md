# Generative Molecular Design Usage Guide

## Overview

Generate novel molecules biased toward desired properties using REINVENT 4 (de novo, scaffold decoration, linker design, molecular optimization). Combine activity, ADMET, drug-likeness, and synthetic accessibility into a multi-objective scoring function. Also covers diffusion-based generation (DiffSMol, DiGress) and JT-VAE alternatives.

## Prerequisites

```bash
pip install reinvent rdkit pytorch
```

REINVENT 4 (Loeffler 2024, AstraZeneca, Apache 2.0).

## Quick Start

Tell the AI agent what to do:
- "Generate 100 de novo molecules optimized for kinase X activity and QED"
- "Decorate this scaffold with 50 novel R-groups, optimize for binding"
- "Design a linker between fragment A and fragment B for PROTAC"
- "Use RL to optimize a lead for activity + low hERG liability"

## Example Prompts

### De novo design with MPO
> "Run REINVENT 4 in de novo mode for 500 steps with this scoring function: 0.4 weight on kinase QSAR, 0.2 QED, 0.2 SA score, 0.2 Tanimoto diversity. Output 100 unique molecules."

### Scaffold decoration
> "Decorate scaffold 'c1ccc(NC(=O)[*:1])cc1' for 200 RL steps. Optimize predicted pIC50 for target X using qsar_model.pkl. Output top 50 with QED > 0.5."

### Lead optimization
> "Optimize lead compound (SMILES given) using REINVENT 4. Maintain Tanimoto >= 0.5 to lead but improve hERG safety. 100 RL steps."

### PROTAC linker design
> "Design 30 linkers between target-side fragment A and E3-side fragment B using REINVENT 4 linker mode. Score by predicted ternary complex stability."

## What the Agent Will Do

1. Set up REINVENT 4 config (TOML) with chosen generator + algorithm.
2. Define scoring function (geometric mean of activity, QED, SA, diversity).
3. Run RL/TL/CL training for N steps.
4. Filter generated molecules to those passing all components.
5. Standardize and deduplicate outputs.
6. Report novelty (Tanimoto to known), drug-likeness, predicted activity.

## Tips

- Always include synthetic accessibility (SA score) to avoid unsynthesizable outputs.
- Use geometric mean (not arithmetic) for combining scoring components.
- Watch for mode collapse: monitor pairwise Tanimoto in generated batch.
- Validate top-N with retrosynthesis (AiZynthFinder).
- PAINS filter as `filter_only=true`; do not reward avoidance.
- For RL: start with 100 steps; increase if score plateaus.

## Related Skills

- chemoinformatics/qsar-modeling - Build scoring models for generation
- chemoinformatics/retrosynthesis - Validate synthetic feasibility
- chemoinformatics/molecular-standardization - Standardize generated SMILES
- chemoinformatics/admet-prediction - ADMET in scoring
- chemoinformatics/substructure-search - PAINS / BRENK filter
- chemoinformatics/scaffold-analysis - Scaffold-aware generation
- chemoinformatics/reaction-enumeration - Combinatorial alternative
