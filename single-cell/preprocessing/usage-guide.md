# Single-Cell Preprocessing - Usage Guide

## Overview

This skill covers quality control, filtering, and normalization for single-cell RNA-seq data using both Seurat (R) and Scanpy (Python). These are essential steps before clustering and downstream analysis.

## Prerequisites

**Python (Scanpy):**
```bash
pip install scanpy matplotlib
```

**R (Seurat):**
```r
install.packages('Seurat')
```

## Quick Start

Ask your AI agent:

> "Run QC on my single-cell data and filter low-quality cells"

> "Normalize my scRNA-seq data and find highly variable genes"

> "Preprocess this 10X data for clustering"

## Example Prompts

### Quality Control
> "Calculate QC metrics including mitochondrial percentage"

> "Show violin plots of QC metrics"

> "What are good filtering thresholds for this dataset?"

### Filtering
> "Filter cells with less than 200 genes or more than 20% mitochondrial"

> "Remove low-quality cells and rarely detected genes"

### Normalization
> "Normalize using log normalization"

> "Run SCTransform on this Seurat object"

> "Normalize to 10,000 counts per cell"

### Feature Selection
> "Find the top 2000 highly variable genes"

> "Show a plot of variable features"

## What the Agent Will Do

1. Calculate QC metrics (gene counts, UMI counts, mito %)
2. Visualize distributions to inform filtering
3. Apply filtering thresholds
4. Normalize and log-transform counts
5. Identify highly variable genes
6. Scale data for PCA

## Tips

- **Store raw counts** before normalization for later use
- **SCTransform** is recommended for Seurat workflows (combines normalize, HVG, scale)
- **Mitochondrial threshold** varies by tissue (5% for PBMCs, 20% for some tissues)
- **Filter doublets first** - high gene counts often indicate doublets
- **Check QC plots** before choosing thresholds - they're dataset-specific
