# Context-Specific Models - Usage Guide

## Overview

Build tissue and condition-specific metabolic models by integrating gene expression data with genome-scale models.

## Prerequisites

```bash
pip install cobra numpy pandas
```

## Quick Start

Tell your AI agent what you want to do:
- "Create a liver-specific model using GTEx expression data"
- "Build a cancer cell model from my RNA-seq data"
- "Integrate my scRNA-seq with the human metabolic model"

## Example Prompts

### Tissue-Specific

> "Create a hepatocyte-specific model from Recon3D"

> "Build tissue models for liver, muscle, and brain"

### Expression Integration

> "Use my TPM data to constrain the E. coli model"

> "Integrate proteomics data with the yeast model"

### Algorithm Choice

> "Apply GIMME to create a context-specific model"

> "Use iMAT to maximize expression-flux agreement"

### Validation

> "Compare fluxes between generic and liver-specific models"

> "Check which pathways are active in my context model"

## What the Agent Will Do

1. Load generic genome-scale model
2. Process expression data (normalize, threshold)
3. Map gene expression to reactions
4. Apply GIMME/iMAT algorithm
5. Validate context model maintains growth
6. Report active/inactive pathways

## Tips

- **Threshold** - 25th percentile is common for "inactive" genes
- **Gene mapping** - Ensure gene IDs match between expression data and model
- **Growth constraint** - Require minimum biomass for biologically meaningful models
- **Validation** - Compare model predictions with known tissue metabolism
- **GPR logic** - Consider AND/OR relationships in gene rules
