---
name: bio-workflows-spatial-pipeline
description: End-to-end spatial transcriptomics workflow for Visium/Xenium data. Covers data loading, preprocessing, spatial analysis, domain detection, and visualization with Squidpy. Use when analyzing spatial transcriptomics data.
tool_type: python
primary_tool: Squidpy
goal_approach_exempt: true
workflow: true
depends_on:
  - spatial-transcriptomics/spatial-data-io
  - spatial-transcriptomics/spatial-preprocessing
  - spatial-transcriptomics/image-analysis
  - spatial-transcriptomics/spatial-deconvolution
  - spatial-transcriptomics/spatial-neighbors
  - spatial-transcriptomics/spatial-statistics
  - spatial-transcriptomics/spatial-domains
  - spatial-transcriptomics/spatial-visualization
qc_checkpoints:
  - platform_fork: "Imaging vs sequencing decided; QC floors and deconvolution-vs-segmentation set accordingly"
  - after_loading: "Spots/cells detected, image aligned"
  - after_qc: "Low-quality spots filtered, genes detected"
  - after_clustering: "Spatial domains correspond to tissue regions"
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, scanpy 1.10+, squidpy 1.3+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Spatial Transcriptomics Pipeline

**"Analyze my spatial transcriptomics data end-to-end"** -> Orchestrate data loading (squidpy/scanpy), QC, normalization, spatial neighbor analysis, spatial statistics, spatial domain detection, and tissue visualization. Composition estimation (deconvolution, spatial-deconvolution) and cell-cell communication (spatial-communication) are deliberately separate steps -- this pipeline hands off to those skills rather than inlining them.

Complete workflow for analyzing Visium, Xenium, or other spatial transcriptomics data.

## The platform-class fork (decide first)

This pipeline branches on platform class before any step. Sequencing/capture data (Visium, Visium HD, Slide-seq, Stereo-seq) are spot/bin MIXTURES of cells: QC on spot counts, normalize knowing that library size partly carries cellularity, then DECONVOLVE composition (spatial-deconvolution) rather than read a spot as one cell. Imaging/in-situ data (Xenium, MERFISH, CosMx) are single molecules: SEGMENT cells first (image-analysis), apply low-count-aware QC floors (an scRNA `min_counts=500` deletes nearly every real imaging cell, whose vector is tens-to-low-hundreds of transcripts), drop on negative-control probe rate, and SKIP deconvolution. The Squidpy+Scanpy path below is written for Visium; the imaging branch is flagged at each step.

## Workflow Overview

```
Spatial data (Space Ranger output)
    |
    v
[1. Load Data] ---------> Read Visium/Xenium
    |
    v
[2. QC & Preprocessing] -> Filter, normalize
    |
    v
[3. Clustering] --------> Standard scRNA-seq clustering
    |
    v
[4. Spatial Analysis] --> Neighbors, statistics
    |
    v
[5. Domain Detection] --> Spatial domains
    |
    v
[6. Visualization] -----> Spatial plots
    |
    v
Annotated spatial data
```

## Primary Path: Squidpy + Scanpy

### Step 1: Load Data

```python
import scanpy as sc
import squidpy as sq
import numpy as np
import matplotlib.pyplot as plt

# Load Visium data (Space Ranger output). squidpy.read provides only visium,
# vizgen, and nanostring -- there is NO sq.read.xenium.
adata = sq.read.visium('spaceranger_output/')

# For Xenium and other imaging platforms use spatialdata_io, which returns a
# SpatialData object preserving the per-transcript molecule table (the cell
# matrix is one table inside it). See spatial-data-io.
# import spatialdata_io as sdio
# sdata = sdio.xenium('xenium_output/')
# adata = sdata.tables['table']  # segmentation-derived cell matrix

print(f'Loaded: {adata.n_obs} spots/cells, {adata.n_vars} genes')
```

### Step 2: Quality Control

```python
# QC metrics. Mito genes are present on Visium but usually OFF-PANEL for imaging
# platforms, so guard the mito calculation rather than assuming MT- genes exist.
has_mito = adata.var_names.str.startswith('MT-').any()
if has_mito:
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
else:
    sc.pp.calculate_qc_metrics(adata, inplace=True)

# Always inspect QC SPATIALLY (a gradient across the section is a technical
# artifact, not biology); violins alone hide it.
sc.pl.spatial(adata, color='total_counts', show=False)
plt.savefig('qc_spatial.pdf')

# Filter. These floors are VISIUM defaults (spot = 1-10-cell mixture) and are
# tissue-dependent. For IMAGING data use low-count-aware floors (~10 transcripts
# per cell, NOT 500) or aggressive filtering deletes nearly every real cell and
# preferentially removes small cells (lymphocytes), biasing composition.
sc.pp.filter_cells(adata, min_counts=500)
sc.pp.filter_genes(adata, min_cells=10)
if has_mito:
    adata = adata[adata.obs.pct_counts_mt < 25, :]

print(f'After QC: {adata.n_obs} spots/cells')
```

### Step 3: Normalization and Clustering

```python
# Store raw counts
adata.layers['counts'] = adata.X.copy()

# Normalize. In spatial data library size partly CARRIES BIOLOGY (Visium total
# counts confound with cells-per-spot and cellularity; imaging total counts with
# cell size), so total-count normalization is a Visium starting point, not a
# universal default -- for imaging consider cell volume/area normalization and
# see spatial-preprocessing before dividing library size out.
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# HVGs
sc.pp.highly_variable_genes(adata, n_top_genes=2000)

# PCA and clustering
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)

# Visualize clusters in space. On Visium these spot clusters are REGIONS/niches,
# NOT cell types -- a spot is a 1-10-cell mixture, so recovering cell-type
# composition needs deconvolution (spatial-deconvolution), not clustering.
sc.pl.spatial(adata, color='leiden', spot_size=1.5)
plt.savefig('clusters_spatial.pdf')
```

### Step 4: Spatial Analysis

```python
# Build spatial neighbors graph. Visium is a hex lattice -> coord_type='grid'
# (n_neighs=6); 'generic' kNN is for imaging point clouds. See spatial-neighbors.
sq.gr.spatial_neighbors(adata, coord_type='grid', n_neighs=6)

# Neighborhood enrichment. The Squidpy permutation null only tests "more adjacent
# than complete spatial randomness" -- two abundant types sharing a compartment
# pass trivially. A positive z is NOT a specific A-B interaction; demand a
# conditional/toroidal null before claiming affinity. See spatial-statistics.
sq.gr.nhood_enrichment(adata, cluster_key='leiden')
sq.pl.nhood_enrichment(adata, cluster_key='leiden')
plt.savefig('nhood_enrichment.pdf')

# Co-occurrence analysis
sq.gr.co_occurrence(adata, cluster_key='leiden')
sq.pl.co_occurrence(adata, cluster_key='leiden')
plt.savefig('co_occurrence.pdf')

# Spatially variable genes. Gate on FDR, not raw I; and a top-Moran gene is
# usually a marker of a spatially-clustered cell TYPE (composition), not a gene
# regulated WITHIN a type -- intersect with non-HVG to find the latter. See
# spatial-statistics.
sq.gr.spatial_autocorr(adata, mode='moran', n_perms=100, n_jobs=4)
moran = adata.uns['moranI']
svg = moran[moran['pval_norm_fdr_bh'] < 0.05].sort_values('I', ascending=False)
print('Spatially autocorrelated genes (FDR<0.05):', svg.head(10).index.tolist())
```

### Step 5: Domain Detection

```python
# Spatial domain detection. Clustering the spatial graph topology ALONE (below)
# is a quick proxy, NOT a real domain method -- it ignores expression and carries
# none of the over-smoothing / spatial-weight-knob / k-as-biological-choice
# framing. For real domains use BANKSY (lambda ~0.8), BayesSpace, or STAGATE and
# tune the spatial weight. See spatial-domains.
sq.gr.spatial_neighbors(adata, coord_type='grid', n_neighs=6)
sc.tl.leiden(adata, resolution=0.3, key_added='spatial_domains',
             adjacency=adata.obsp['spatial_connectivities'])

# Visualize domains
sc.pl.spatial(adata, color='spatial_domains', spot_size=1.5)
plt.savefig('spatial_domains.pdf')

# Compare transcriptomic vs spatial clusters
sc.pl.spatial(adata, color=['leiden', 'spatial_domains'], ncols=2)
plt.savefig('clusters_comparison.pdf')
```

### Step 6: Visualization

```python
# Gene expression in space
genes = ['EPCAM', 'VIM', 'PTPRC', 'COL1A1']
sc.pl.spatial(adata, color=genes, ncols=2, spot_size=1.5, cmap='viridis')
plt.savefig('marker_genes_spatial.pdf')

# Cluster markers in space. On Visium these are markers of spot REGIONS (mixtures),
# not of pure cell types; for cell-type-level signal deconvolve first.
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5)
plt.savefig('cluster_markers.pdf')

# Save
adata.write('spatial_analyzed.h5ad')
```

## Complete Workflow Script

```python
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import os

# Configuration
data_dir = 'spaceranger_output'
output_dir = 'spatial_results'
os.makedirs(output_dir, exist_ok=True)
os.makedirs(f'{output_dir}/plots', exist_ok=True)

# Load
print('Loading data...')
adata = sq.read.visium(data_dir)
print(f'Loaded: {adata.n_obs} spots, {adata.n_vars} genes')

# QC (Visium defaults; for imaging use low-count-aware floors and skip mito)
print('QC filtering...')
has_mito = adata.var_names.str.startswith('MT-').any()
if has_mito:
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
else:
    sc.pp.calculate_qc_metrics(adata, inplace=True)
sc.pp.filter_cells(adata, min_counts=500)
sc.pp.filter_genes(adata, min_cells=10)
if has_mito:
    adata = adata[adata.obs.pct_counts_mt < 25, :]
print(f'After QC: {adata.n_obs} spots')

# Normalize and cluster
print('Processing...')
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.leiden(adata, resolution=0.5)

# Spatial analysis (Visium hex -> coord_type='grid'; nhood z and top-Moran genes
# need the caveats from Step 4 before interpretation)
print('Spatial analysis...')
sq.gr.spatial_neighbors(adata, coord_type='grid', n_neighs=6)
sq.gr.nhood_enrichment(adata, cluster_key='leiden')
sq.gr.spatial_autocorr(adata, mode='moran', n_perms=100)

# Plots
print('Creating plots...')
sc.pl.spatial(adata, color='leiden', spot_size=1.5, save='_clusters.pdf')
sq.pl.nhood_enrichment(adata, cluster_key='leiden', save='_nhood.pdf')

# Save
adata.write(f'{output_dir}/spatial_analyzed.h5ad')
print(f'Results saved to {output_dir}/')
```

## Related Skills

- spatial-transcriptomics/spatial-data-io - Loading formats
- spatial-transcriptomics/spatial-preprocessing - QC details
- spatial-transcriptomics/spatial-statistics - Moran's I, co-occurrence
- spatial-transcriptomics/spatial-domains - Domain detection methods
- spatial-transcriptomics/spatial-deconvolution - Cell type estimation
