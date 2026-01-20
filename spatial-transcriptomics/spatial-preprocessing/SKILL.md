---
name: bio-spatial-transcriptomics-spatial-preprocessing
description: Quality control, filtering, normalization, and feature selection for spatial transcriptomics data. Calculate QC metrics, filter spots/cells, normalize counts, and identify highly variable genes.
tool_type: python
primary_tool: squidpy
---

# Spatial Preprocessing

QC, filtering, normalization, and feature selection for spatial data.

## Required Imports

```python
import squidpy as sq
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
```

## Calculate QC Metrics

```python
# Calculate standard QC metrics
sc.pp.calculate_qc_metrics(adata, inplace=True)

# View QC columns
print(adata.obs[['total_counts', 'n_genes_by_counts']].describe())
print(adata.var[['total_counts', 'n_cells_by_counts']].describe())
```

## Calculate Mitochondrial Content

```python
# Mark mitochondrial genes
adata.var['mt'] = adata.var_names.str.startswith('MT-')

# Calculate percent mitochondrial
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
print(f"Mean MT%: {adata.obs['pct_counts_mt'].mean():.1f}")
```

## Visualize QC Metrics on Tissue

```python
# Plot QC metrics spatially
sq.pl.spatial_scatter(adata, color=['total_counts', 'n_genes_by_counts', 'pct_counts_mt'], ncols=3)

# Or with Scanpy
sc.pl.spatial(adata, color=['total_counts', 'n_genes_by_counts'], spot_size=1.5)
```

## QC Metric Distributions

```python
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
axes[0].hist(adata.obs['total_counts'], bins=50)
axes[0].set_xlabel('Total counts')
axes[1].hist(adata.obs['n_genes_by_counts'], bins=50)
axes[1].set_xlabel('Genes detected')
axes[2].hist(adata.obs['pct_counts_mt'], bins=50)
axes[2].set_xlabel('MT %')
plt.tight_layout()
```

## Filter Spots

```python
# Filter based on QC metrics
print(f'Before filtering: {adata.n_obs} spots')

# Minimum counts and genes
sc.pp.filter_cells(adata, min_counts=500)
sc.pp.filter_cells(adata, min_genes=200)

# Maximum mitochondrial content
adata = adata[adata.obs['pct_counts_mt'] < 20].copy()

print(f'After filtering: {adata.n_obs} spots')
```

## Filter Genes

```python
# Remove genes detected in few spots
print(f'Before filtering: {adata.n_vars} genes')
sc.pp.filter_genes(adata, min_cells=10)
print(f'After filtering: {adata.n_vars} genes')
```

## Normalization

```python
# Store raw counts
adata.layers['counts'] = adata.X.copy()

# Normalize to median total counts
sc.pp.normalize_total(adata, target_sum=1e4)

# Log transform
sc.pp.log1p(adata)
```

## SCTransform-like Normalization

```python
# Pearson residuals normalization (similar to SCTransform)
# Requires raw counts
adata_raw = adata.copy()
adata_raw.X = adata_raw.layers['counts']

sc.experimental.pp.normalize_pearson_residuals(adata_raw)
adata.layers['pearson'] = adata_raw.X.copy()
```

## Highly Variable Genes

```python
# Find HVGs
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3', layer='counts')

# View HVG stats
print(f"Found {adata.var['highly_variable'].sum()} HVGs")
sc.pl.highly_variable_genes(adata)
```

## Spatially Variable Genes

```python
# Compute spatial neighbors first
sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=6)

# Find spatially variable genes using Moran's I
sq.gr.spatial_autocorr(adata, mode='moran', genes=adata.var_names[:1000])

# Get top spatially variable genes
svg = adata.uns['moranI'].sort_values('I', ascending=False)
print('Top spatially variable genes:')
print(svg.head(20))
```

## Combine HVG and SVG

```python
# Get union of highly variable and spatially variable genes
hvg = set(adata.var_names[adata.var['highly_variable']])
svg_top = set(adata.uns['moranI'].head(500).index)
selected_genes = hvg | svg_top

print(f'HVG: {len(hvg)}, SVG: {len(svg_top)}, Union: {len(selected_genes)}')

# Subset to selected genes for downstream
adata_subset = adata[:, list(selected_genes)].copy()
```

## Scale Data

```python
# Scale for PCA (use log-normalized data)
sc.pp.scale(adata, max_value=10)
```

## PCA

```python
# Run PCA
sc.tl.pca(adata, n_comps=50)

# Variance explained
sc.pl.pca_variance_ratio(adata, n_pcs=50)
```

## Complete Preprocessing Pipeline

```python
import squidpy as sq
import scanpy as sc

# Load data
adata = sq.read.visium('spaceranger_output/')

# QC
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# Filter
sc.pp.filter_cells(adata, min_counts=1000)
sc.pp.filter_cells(adata, min_genes=500)
adata = adata[adata.obs['pct_counts_mt'] < 20].copy()
sc.pp.filter_genes(adata, min_cells=10)

# Normalize
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# HVGs
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3', layer='counts')

# Scale and PCA
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50)

print(f'Preprocessed: {adata.n_obs} spots, {adata.n_vars} genes')
adata.write_h5ad('preprocessed.h5ad')
```

## Related Skills

- **spatial-data-io** - Load spatial data
- **spatial-neighbors** - Build spatial graphs
- **single-cell/preprocessing** - Non-spatial preprocessing
