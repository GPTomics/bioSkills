# Single-Cell Data I/O - Usage Guide

## Overview

This skill covers reading, writing, creating, and converting single-cell objects across AnnData (Python), Seurat (R), and SingleCellExperiment (R). It emphasizes the decisions that prevent silent data loss: keeping the raw 10X matrix, choosing stable gene identifiers, picking a storage format, and converting between Python and R without dropping layers or transposing the matrix incorrectly.

## Prerequisites

**Python (Scanpy/AnnData):**
```bash
pip install scanpy anndata
# multimodal: pip install muon mudata
```

**R (Seurat + conversion):**
```r
install.packages('Seurat')
# Conversion (prefer maintained tools; SeuratDisk is abandoned):
remotes::install_github('scverse/anndataR')          # pure-R h5ad/zarr I/O + conversion
BiocManager::install('zellkonverter')                # SCE <-> AnnData
remotes::install_github('cellgeni/schard')           # robust pure-R h5ad reading
```

## Quick Start

Tell the AI agent what is needed:
- "Load the raw 10X matrix and keep the antibody-capture features"
- "Create an AnnData object from this count matrix and store raw counts"
- "Convert this h5ad to a Seurat object without losing the UMAP and layers"

## Example Prompts

### Loading Data
> "Read the raw_feature_bc_matrix folder using Ensembl gene IDs"

> "Load this Cell Ranger h5 and report cells x genes"

> "Load the 10X output but keep CRISPR guide and antibody features"

### Creating Objects
> "Build an AnnData from this matrix, put integer counts in a counts layer"

> "Create a Seurat v5 object from this sparse matrix with min.cells 3"

### Converting
> "Convert this AnnData to a SingleCellExperiment, keeping reducedDims and raw"

> "Move this Seurat object to h5ad for Python and verify no layers were dropped"

> "Why did my gene and cell axes swap after conversion?"

## What the Agent Will Do

1. Identify the input format and whether the raw (unfiltered) 10X matrix is available
2. Choose gene identifiers (Ensembl IDs for reproducibility) and whether to retain non-GEX features
3. Read into AnnData or Seurat, placing counts/normalized/embeddings in conventional slots
4. For conversion, select a maintained converter, apply the transpose, and remap metadata axes
5. Diff slot inventories before and after conversion to confirm nothing was silently dropped
6. Write to the appropriate format (h5ad/zarr/RDS/h5mu) and keep the original

## Tips

- **Keep the raw matrix** - EmptyDrops, SoupX, CellBender, and DecontX all need the unfiltered Cell Ranger output; filtered-only storage is irreversible.
- **Filtered is not decontaminated** - Cell Ranger filtered output is cell-CALLED, not ambient-corrected; that is a separate step.
- **Use gene IDs for joins** - gene symbols are non-unique and change across annotation releases; Ensembl IDs are stable.
- **Mind the transpose** - AnnData is cells x genes; Seurat/SCE are genes x cells. Conversion transposes the matrix AND swaps which axis the metadata annotates.
- **Conversion is lossy by default** - layers, obsp/varp, nested uns, and categoricals can vanish; diff slots before and after.
- **Avoid SeuratDisk** - abandoned since 2023 and broken on Seurat v5; prefer anndataR, zellkonverter, or schard.
- **Keep matrices sparse** - dense materialization of a large object exhausts memory; check `scipy.sparse.issparse(adata.X)`.
- **Set drop_single_values=FALSE** - sceasy deletes constant metadata columns by default, losing single-sample batch labels.

## Related Skills

- single-cell/preprocessing - QC, normalization, and HVG selection after loading
- single-cell/doublet-detection - per-sample doublet calling on raw counts after loading
- single-cell/clustering - dimensionality reduction and clustering on the loaded object
- single-cell/multimodal-integration - MuData/h5mu handling for CITE-seq and Multiome
- spatial-transcriptomics/spatial-data-io - SpatialData/zarr I/O for spatial omics
- workflows/scrnaseq-pipeline - end-to-end scRNA-seq pipeline that starts from data loading
