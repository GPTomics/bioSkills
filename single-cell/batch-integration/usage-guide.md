# Batch Integration Usage Guide

Integrate multiple scRNA-seq datasets to remove batch effects while preserving biological variation.

## Prerequisites

```bash
# Python
pip install scanpy harmonypy scvi-tools

# R
install.packages('Seurat')
install.packages('harmony')
BiocManager::install('batchelor')
```

## Quick Start

### Harmony (R) - Fast, Most Common

```r
library(Seurat)
library(harmony)

# Merge datasets
merged <- merge(sample1, y = list(sample2, sample3), add.cell.ids = c('S1', 'S2', 'S3'))

# Standard preprocessing
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged)
merged <- RunPCA(merged)

# Run Harmony
merged <- RunHarmony(merged, group.by.vars = 'orig.ident')

# Continue with integrated embeddings
merged <- RunUMAP(merged, reduction = 'harmony', dims = 1:30)
merged <- FindNeighbors(merged, reduction = 'harmony', dims = 1:30)
merged <- FindClusters(merged, resolution = 0.5)
```

### Harmony (Python)

```python
import scanpy as sc
import harmonypy

adata = sc.read_h5ad('merged.h5ad')
sc.pp.pca(adata)

# Run Harmony
adata.obsm['X_pca_harmony'] = harmonypy.run_harmony(
    adata.obsm['X_pca'], adata.obs, 'batch'
).Z_corr.T

sc.pp.neighbors(adata, use_rep='X_pca_harmony')
sc.tl.umap(adata)
```

## Method Comparison

| Method | Speed | Best For |
|--------|-------|----------|
| Harmony | Fast | Most use cases |
| scVI | Moderate | Large datasets, deep learning |
| Seurat CCA | Moderate | Conserved biology |
| fastMNN | Fast | MNN-based correction |

## scVI (Deep Learning)

```python
import scvi

scvi.model.SCVI.setup_anndata(adata, batch_key='batch')
model = scvi.model.SCVI(adata)
model.train()

adata.obsm['X_scVI'] = model.get_latent_representation()
sc.pp.neighbors(adata, use_rep='X_scVI')
```

## Seurat CCA/RPCA

```r
# Split by batch
samples <- SplitObject(merged, split.by = 'batch')

# Find anchors
anchors <- FindIntegrationAnchors(object.list = samples, dims = 1:30)

# Integrate
integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

# Use 'integrated' assay for downstream
DefaultAssay(integrated) <- 'integrated'
```

## Evaluating Integration

### Visual Assessment

- UMAP should show mixing of batches within cell types
- Cell types should cluster together across batches

### Quantitative Metrics

```python
# kBET - batch mixing within neighborhoods
# LISI - local inverse Simpson index
# Silhouette - cluster separation
```

## Tips

- Always preprocess each batch separately before merging
- Check that cell types are represented across batches
- Use batch as covariate in DE analysis, not integrated values
- For DE, use original counts, not batch-corrected values

## See Also

- [Harmony paper](https://www.nature.com/articles/s41592-019-0619-0)
- [scVI documentation](https://docs.scvi-tools.org/)
