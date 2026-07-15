# Generative Molecular Design Usage Guide

## Overview

Generate novel molecules biased toward desired properties using REINVENT 4 (de novo, scaffold decoration, linker design, molecular optimization). Combine activity, ADMET, drug-likeness, and synthetic accessibility into a multi-objective scoring function. Also covers diffusion-based generation (DiffSMol, DiGress) and JT-VAE alternatives.

## Prerequisites

```bash
git clone https://github.com/MolecularAI/REINVENT4.git
cd REINVENT4
pip install -e .
```

Follow the selected REINVENT 4 release's installation instructions and use its matching example configs. `pip install reinvent` does not install REINVENT 4.

## Quick Start

Tell the AI agent what to do:
- "Generate 100 de novo molecules optimized for kinase X activity and QED"
- "Decorate this scaffold with 50 novel R-groups, optimize for binding"
- "Design a linker between fragment A and fragment B for PROTAC"
- "Use RL to optimize a lead for activity + low hERG liability"

## Example Prompts

### De novo design with MPO
> "Run a release-matched REINVENT 4 de novo pilot. Justify each scoring component from project data, monitor learning/diversity curves, and post-process the stage CSV."

### Scaffold decoration
> "Decorate scaffold 'c1ccc(NC(=O)[*:1])cc1' for 200 RL steps. Optimize predicted pIC50 for target X using qsar_model.pkl. Output top 50 with QED > 0.5."

### Lead optimization
> "Optimize lead compound (SMILES given) using REINVENT 4. Maintain Tanimoto >= 0.5 to lead but improve hERG safety. 100 RL steps."

### PROTAC linker design
> "Design 30 linkers between target-side fragment A and E3-side fragment B using REINVENT 4 linker mode. Score by predicted ternary complex stability."

## What the Agent Will Do

1. Set up REINVENT 4 config (TOML) with chosen generator + algorithm.
2. Define and validate a project-specific scoring function; aggregation may be geometric, arithmetic, or otherwise justified.
3. Run RL/TL/CL training for N steps.
4. Annotate and prioritize generated molecules using project-defined decision rules.
5. Standardize and deduplicate outputs.
6. Report novelty (Tanimoto to known), drug-likeness, predicted activity.

## Tips

- Treat SA score as a heuristic annotation, not proof of a route; validate selected molecules with reaction- or route-based methods.
- Choose geometric versus arithmetic aggregation from the desired compensation behavior and validation results.
- REINVENT writes a live CSV per stage and the configured checkpoint at stage termination; it does not emit per-iteration `.smi` and checkpoint files by default.
- Watch for mode collapse: monitor pairwise Tanimoto in generated batch.
- Validate top-N with retrosynthesis (AiZynthFinder).
- REINVENT 4 `CustomAlerts` zeroes matches before score aggregation; there is no `filter_only` setting. Treat PAINS as triage alerts unless the project defines a justified exclusion rule.
- Treat 100 RL steps, the prompt weights, and prompt similarity cutoffs as repository starting values only; tune them against learning curves, diversity, oracle validity, and project objectives.

## Related Skills

- chemoinformatics/qsar-modeling - Build scoring models for generation
- chemoinformatics/retrosynthesis - Validate synthetic feasibility
- chemoinformatics/molecular-standardization - Standardize generated SMILES
- chemoinformatics/admet-prediction - ADMET in scoring
- chemoinformatics/substructure-search - PAINS / BRENK filter
- chemoinformatics/scaffold-analysis - Scaffold-aware generation
- chemoinformatics/reaction-enumeration - Combinatorial alternative
