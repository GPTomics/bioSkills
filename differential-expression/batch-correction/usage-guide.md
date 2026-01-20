# Batch Correction - Usage Guide

## Overview

Batch effects are technical variations between experimental batches that can confound biological signals. Different correction methods suit different analysis goals.

## Prerequisites

```r
BiocManager::install(c('sva', 'limma', 'DESeq2'))
install.packages('harmony')  # for single-cell
```

## Quick Start

Tell your AI agent what you want to do:
- "Remove batch effects from my RNA-seq data using ComBat-Seq"
- "Add batch as a covariate in DESeq2 analysis"
- "Estimate surrogate variables for unknown batch effects"

## Example Prompts

### Design-Based Correction
> "Include batch in my DESeq2 design formula"

> "Model batch effects as a covariate in edgeR"

### Count Correction
> "Apply ComBat-Seq to my raw count matrix"

> "Correct counts for batch before clustering"

### Unknown Batches
> "Use SVA to find hidden batch effects in my data"

> "Estimate surrogate variables and include in DE analysis"

### Visualization
> "Make a PCA plot colored by batch before and after correction"

> "Assess batch effect severity in my samples"

## What the Agent Will Do

1. Visualize batch effects with PCA colored by batch
2. Choose appropriate correction method based on downstream analysis
3. Apply correction (design formula for DE, ComBat-Seq for visualization)
4. Validate correction with post-correction PCA

## Method Selection

| Method | Input | Use For |
|--------|-------|---------|
| DESeq2 design formula | Raw counts | DE analysis (preferred) |
| ComBat-Seq | Raw counts | Visualization, clustering |
| ComBat | Normalized | Visualization, ML |
| limma removeBatchEffect | Normalized | Visualization only |
| SVA | Normalized | Unknown batch sources |
| Harmony | Embeddings | Single-cell integration |

## Tips

- For differential expression: Include batch in the design formula, don't correct counts
- For visualization/clustering: Use corrected values
- Known vs unknown batches: Use SVA when batch sources are unknown
- Confounding: Batch perfectly correlated with condition is unfixable
- Over-correction: Can remove biological signal if not careful
- Balanced design: Best results when conditions are spread across batches
