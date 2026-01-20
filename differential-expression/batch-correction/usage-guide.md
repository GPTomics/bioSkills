# Batch Correction Usage Guide

## Overview

Batch effects are technical variations between experimental batches that can confound biological signals. Different correction methods suit different analysis goals.

## Method Selection

| Method | Input | Use For |
|--------|-------|---------|
| DESeq2 design formula | Raw counts | DE analysis (preferred) |
| ComBat-Seq | Raw counts | Visualization, clustering |
| ComBat | Normalized | Visualization, ML |
| limma removeBatchEffect | Normalized | Visualization only |
| SVA | Normalized | Unknown batch sources |
| Harmony | Embeddings | Single-cell integration |

## Quick Start Prompts

- "Remove batch effects from my RNA-seq data using ComBat-Seq"
- "Add batch as a covariate in DESeq2 analysis"
- "Estimate surrogate variables for unknown batch effects"
- "Visualize batch effects with PCA before and after correction"

## Key Principle

**For differential expression**: Include batch in the design formula, don't correct counts
**For visualization/clustering**: Use corrected values

## Workflow

1. **Visualize** - PCA colored by batch to assess severity
2. **Choose method** - Based on downstream analysis
3. **Correct/model** - Apply appropriate method
4. **Validate** - PCA should show condition, not batch, as main driver

## Requirements

```r
BiocManager::install(c('sva', 'limma', 'DESeq2'))
install.packages('harmony')  # for single-cell
```

## Key Considerations

- **Known vs unknown batches** - SVA for unknown
- **Confounding** - Batch perfectly correlated with condition is unfixable
- **Over-correction** - Can remove biological signal
- **Balanced design** - Best when conditions spread across batches

## Related Skills

- **differential-expression/deseq2-basics** - DE analysis
- **single-cell/clustering-annotation** - SC integration
