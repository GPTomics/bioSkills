# Single-Cell Data I/O - Usage Guide

## Overview

This skill covers reading, writing, and creating single-cell data objects using both Seurat (R) and Scanpy (Python). Use it for loading 10X Genomics data, importing/exporting files, and managing cell and gene metadata.

## Prerequisites

**Python (Scanpy):**
```bash
pip install scanpy anndata
```

**R (Seurat):**
```r
install.packages('Seurat')
# For format conversion:
remotes::install_github('mojaveazure/seurat-disk')
```

## Quick Start

Ask your AI agent:

> "Load my 10X data into Scanpy"

> "Create a Seurat object from this count matrix"

> "Convert my h5ad file to Seurat format"

## Example Prompts

### Loading Data
> "Read the 10X filtered_feature_bc_matrix folder"

> "Load this h5ad file and show the cell count"

> "Import the cellranger h5 output"

### Creating Objects
> "Create an AnnData object from this count matrix CSV"

> "Make a Seurat object from this sparse matrix"

### Metadata
> "Add sample labels to the cell metadata"

> "Show me all the cell metadata columns"

### Saving/Converting
> "Save this AnnData as h5ad"

> "Convert this Seurat object to h5ad for use in Python"

## What the Agent Will Do

1. Identify data format and choose appropriate reader
2. Load data into Seurat object or AnnData object
3. Add requested metadata
4. Save in requested format

## Tool Comparison

| Task | Scanpy (Python) | Seurat (R) |
|------|-----------------|------------|
| Read 10X | `sc.read_10x_mtx()` | `Read10X()` |
| Create object | `ad.AnnData()` | `CreateSeuratObject()` |
| Native format | `.h5ad` | `.rds` |
| Access counts | `adata.X` | `LayerData(obj, layer='counts')` |

## Tips

- **Use h5ad for Python workflows** - Native AnnData format, efficient for large datasets
- **Use RDS for R workflows** - Preserves all Seurat object structure
- **Store raw counts** - Use `adata.raw` or `adata.layers['counts']` before normalization
- **Seurat v5 uses layers** - Not slots; use `LayerData()` instead of `GetAssayData()`
- **SeuratDisk for conversion** - Required for Seurat <-> AnnData conversion
