# Single-Cell RNA-seq Pipeline - Usage Guide

## Overview

This workflow processes single-cell RNA-seq data from 10X Genomics Cell Ranger output to annotated cell types. It supports both Seurat (R) and Scanpy (Python) implementations.

## Prerequisites

```r
# R packages
install.packages(c('Seurat', 'ggplot2', 'dplyr'))
BiocManager::install(c('scDblFinder', 'SingleCellExperiment'))
```

```bash
# Python packages
pip install scanpy scrublet anndata
```

## Quick Start

Tell your AI agent what you want to do:
- "Analyze my 10X single-cell data from start to finish"
- "Run the scRNA-seq pipeline on my PBMC data"
- "Cluster my single-cell data and find marker genes"

## Example Prompts

### Starting from Cell Ranger output
> "I have filtered_feature_bc_matrix from Cell Ranger, analyze it"

> "Load my 10X data and perform QC filtering"

> "Process my scRNA-seq with Scanpy instead of Seurat"

### Analysis steps
> "Remove doublets from my single-cell data"

> "Find clusters at different resolutions"

> "What are the marker genes for cluster 3?"

### Annotation
> "Help me annotate cell types based on marker genes"

> "Use SingleR for automated cell type annotation"

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| 10X data | filtered_feature_bc_matrix/ | Directory with matrix.mtx, barcodes.tsv, features.tsv |
| 10X H5 | .h5 | Alternatively, the HDF5 format |

## What the Workflow Does

1. **Load Data** - Read 10X matrix into Seurat/AnnData
2. **QC Filtering** - Remove low-quality cells and doublets
3. **Normalization** - SCTransform or log-normalization
4. **Feature Selection** - Find highly variable genes
5. **Dimension Reduction** - PCA and UMAP
6. **Clustering** - Graph-based clustering
7. **Markers** - Find cluster-specific genes
8. **Annotation** - Assign cell type labels

## Seurat vs Scanpy

| Feature | Seurat | Scanpy |
|---------|--------|--------|
| Language | R | Python |
| Speed | Fast | Faster for large data |
| Memory | Moderate | Lower |
| Ecosystem | Bioconductor | Python ML stack |
| Best for | General analysis | Large datasets, integration |

## Tips

- **Cell numbers**: Expect 1,000-20,000 cells from typical 10X run
- **Genes per cell**: 200-5,000 is typical; very high may be doublets
- **Mitochondrial**: >20% suggests dying cells
- **Resolution**: Start at 0.5, adjust based on cluster quality
- **Annotation**: Check canonical markers for your tissue type
