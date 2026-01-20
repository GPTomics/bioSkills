# Batch Integration - Usage Guide

## Overview

Batch integration removes technical variation between samples, experiments, or technologies while preserving biological differences, enabling meaningful comparison across datasets.

## Prerequisites

```bash
# Python
pip install scanpy harmonypy scvi-tools
```

```r
# R
install.packages('Seurat')
install.packages('harmony')
BiocManager::install('batchelor')
```

## Quick Start

Tell your AI agent what you want to do:
- "Integrate my samples to remove batch effects"
- "Run Harmony on my merged Seurat object"
- "Combine datasets from different experiments"

## Example Prompts

### Integration
> "Merge my samples and run Harmony integration"
> "Use scVI to integrate these batches"
> "Run Seurat CCA integration on my samples"

### Assessment
> "Show batch mixing on the UMAP"
> "Calculate integration metrics (LISI, kBET)"
> "Are the batches well mixed within cell types?"

### Comparison
> "Compare Harmony vs scVI integration"
> "Which method preserves cell type separation best?"
> "Test different integration parameters"

### Downstream
> "Cluster the integrated data"
> "Find markers using the integrated representation"
> "Run differential expression between conditions"

## What the Agent Will Do

1. Merge datasets with batch labels
2. Preprocess each batch appropriately
3. Run chosen integration method
4. Generate integrated low-dimensional representation
5. Visualize batch mixing
6. Cluster using integrated embeddings
7. Assess integration quality

## Method Comparison

| Method | Speed | Best For |
|--------|-------|----------|
| Harmony | Fast | Most use cases |
| scVI | Moderate | Large datasets, deep learning |
| Seurat CCA | Moderate | Conserved biology |
| fastMNN | Fast | MNN-based correction |

## Evaluating Integration

### Visual Assessment
- UMAP should show mixing of batches within cell types
- Cell types should cluster together across batches

### Quantitative Metrics
- kBET: batch mixing within neighborhoods
- LISI: local inverse Simpson index
- Silhouette: cluster separation

## Tips

- **Preprocess each batch separately** before merging
- **Check cell type representation** - ensure types are present across batches
- **Use batch as covariate in DE** - not integrated values
- **Keep original counts for DE** - use raw counts, not batch-corrected
- **Validate integration** - cell types should mix, not batch artifacts
- **Harmony is a good default** - fast and works well for most cases
