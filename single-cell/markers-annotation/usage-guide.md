# Marker Genes and Cell Type Annotation - Usage Guide

## Overview

This skill covers finding differentially expressed marker genes and annotating cell types in single-cell RNA-seq data using both Seurat (R) and Scanpy (Python).

## Prerequisites

**Python (Scanpy):**
```bash
pip install scanpy pandas matplotlib
```

**R (Seurat):**
```r
install.packages(c('Seurat', 'dplyr'))
# Optional for better DE:
BiocManager::install('MAST')
```

## Quick Start

Ask your AI agent:

> "Find marker genes for each cluster"

> "What cell types are in my data?"

> "Show a dot plot of marker expression"

## Example Prompts

### Finding Markers
> "Find differentially expressed genes for cluster 0"

> "What genes distinguish cluster 1 from cluster 2?"

> "Show top 10 markers for each cluster"

### Visualization
> "Create a dot plot of these marker genes"

> "Show a heatmap of the top markers"

> "Plot CD3D expression on UMAP"

### Annotation
> "Annotate clusters based on these markers"

> "Score cells for T cell signature genes"

> "Label the clusters with cell type names"

### Export
> "Export all markers to CSV"

> "Save the top 20 markers per cluster"

## What the Agent Will Do

1. Run differential expression test
2. Filter for significant markers
3. Visualize marker expression patterns
4. Suggest cell type annotations
5. Add annotations to object
6. Export results

## Tips

- **Use wilcoxon** for quick marker detection
- **Use MAST** (Seurat) for rigorous DE analysis
- **Check known markers** before automatic annotation
- **Lower min.pct/logfc thresholds** if few markers found
- **Visualize before annotating** - use dot plots and UMAPs
- **Store raw counts** - needed for some DE methods
