# Metabolic Modeling Pipeline - Usage Guide

## Overview

This workflow builds genome-scale metabolic models (GEMs) from genomic data and uses them to predict cellular phenotypes. It covers automated reconstruction, iterative curation, flux balance analysis, and applications like gene essentiality prediction and context-specific modeling.

## Prerequisites

```bash
pip install cobra carveme memote escher pandas numpy matplotlib

conda install -c bioconda diamond
```

**Optional solvers (faster for large models):**
```bash
conda install -c gurobi gurobi
conda install -c ibmdecisionoptimization cplex
```

## Quick Start

Tell your AI agent what you want to do:
- "Build a metabolic model from my genome annotation"
- "Run FBA to predict growth rate on glucose minimal medium"
- "Find essential genes in my metabolic model"
- "Gap-fill my model to grow on M9 medium"

## Example Prompts

### Model building
> "Create a genome-scale metabolic model from my E. coli protein FASTA"

> "Reconstruct a GEM for this Pseudomonas genome"

### Model validation
> "Run memote QC on my metabolic model"

> "Check for blocked reactions and dead-end metabolites"

### Flux analysis
> "Predict growth rate and flux distribution on glucose"

> "Run FVA to find flux ranges for glycolysis reactions"

> "What's the maximum theoretical ethanol yield?"

### Gene essentiality
> "Predict essential genes in my model"

> "Find synthetic lethal gene pairs"

### Context-specific
> "Build a liver-specific model using my RNA-seq data"

> "Constrain the model to match measured uptake rates"

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| Protein sequences | FASTA | Annotated proteome for reconstruction |
| Existing model | SBML/JSON | For analysis or curation |
| Media composition | Dict/TSV | Exchange reaction bounds |
| Expression data (optional) | TSV | Gene-level TPM for context models |

## What the Agent Will Do

1. **Reconstruction** - Generates draft model from protein sequences using CarveMe
2. **Validation** - Runs memote QC to identify model issues
3. **Curation** - Gap-fills for growth, fixes dead-ends, adds missing GPRs
4. **FBA Analysis** - Predicts optimal growth and flux distribution
5. **Applications** - Gene essentiality, context-specific models, yield prediction

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| Memote score threshold | 50% | Minimum for usable model |
| Growth threshold | 0.01 h^-1 | Minimum viable growth |
| FVA optimality fraction | 0.9 | Allow 90% of max growth |
| Essentiality threshold | 10% WT | Below = essential |
| Expression percentile | 25th | Context model cutoff |

## Model Quality Checklist

A well-curated model should have:
- [ ] Growth rate in realistic range (0.1-1.0 h^-1 for bacteria)
- [ ] Memote score >50% (ideally >70%)
- [ ] <100 blocked reactions
- [ ] <50 dead-end metabolites
- [ ] Essential genes overlap >70% with experimental data
- [ ] Key pathways (glycolysis, TCA, etc.) carry flux

## Tips

- **Start with CarveMe**: Fastest and most automated reconstruction
- **Gap-fill iteratively**: Fix one issue at a time, re-test growth
- **Validate against data**: Compare predictions to experimental phenotypes
- **Use commercial solvers**: Gurobi/CPLEX are much faster than GLPK
- **Document changes**: Keep track of manual curation steps
- **Constrain realistically**: Set uptake bounds based on experimental data
- **Test edge cases**: Ensure model behaves correctly under different conditions
