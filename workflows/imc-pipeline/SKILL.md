---
name: bio-workflows-imc-pipeline
description: End-to-end imaging mass cytometry workflow from raw acquisitions to spatial cell analysis. Orchestrates image preprocessing, segmentation, phenotyping, and spatial statistics. Use when analyzing imaging mass cytometry data end-to-end.
tool_type: python
primary_tool: steinbock
workflow: true
depends_on:
  - imaging-mass-cytometry/data-preprocessing
  - imaging-mass-cytometry/cell-segmentation
  - imaging-mass-cytometry/phenotyping
  - imaging-mass-cytometry/spatial-analysis
  - imaging-mass-cytometry/differential-analysis
  - imaging-mass-cytometry/interactive-annotation
  - imaging-mass-cytometry/quality-metrics
---

## Version Compatibility

Reference examples tested with: Cellpose 3.0+, anndata 0.10+, matplotlib 3.8+, numpy 1.26+, pandas 2.2+, scanpy 1.10+, scvi-tools 1.1+, squidpy 1.3+, steinbock 0.16+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Imaging Mass Cytometry Pipeline

**"Process my imaging mass cytometry data from images to spatial analysis"** -> Orchestrate image preprocessing (steinbock), cell segmentation (Cellpose), phenotyping (FlowSOM/scanpy), spatial neighborhood analysis (squidpy), and tissue community detection.

## Pipeline Overview

```
Raw MCD/TIFF Files ──> Image Processing ──> Cell Masks
                                                 │
                                                 ▼
                ┌─────────────────────────────────────────────┐
                │              imc-pipeline                   │
                ├─────────────────────────────────────────────┤
                │  1. Data Preprocessing (spillover, hot px)  │
                │  2. Cell Segmentation (Cellpose/Mesmer)     │
                │  3. Single-cell Quantification              │
                │  4. Clustering & Phenotyping                │
                │  5. Spatial Analysis                        │
                │  6. Visualization                           │
                └─────────────────────────────────────────────┘
                                                 │
                                                 ▼
                    Cell Types + Spatial Neighborhoods
```

## Decisions Threaded Through This Pipeline

Four reframes govern every stage and are detailed in the depended-on skills: IMC pixels are integer ion COUNTS (arcsinh cofactor 1, not the suspension-CyTOF 5), and spillover is spatial so it must be NNLS-compensated before segmentation; segmentation is the largest irreversible error source, so impossible double-positives are a QC alarm, not biology; a spatial interaction is a hypothesis test whose null silently decides whether the result is real or a density artifact; and the experimental unit is the patient, not the cell, so cross-condition tests aggregate to patients before testing.

## Complete steinbock Workflow

### Step 1: Setup and Preprocessing

```bash
# generate the panel template; edit the keep column before extracting
steinbock preprocess imc panel

# extract per-channel TIFFs (keep-filtered, panel-ordered) with hot-pixel removal
# (--hpf is a signed 8-neighbor difference; 50 is a count, tune to dynamic range)
steinbock preprocess imc images --hpf 50

# channel spillover is compensated with NNLS (CATALYST/cytomapper, R) on the pixel images
# BEFORE segmentation when spatial analysis is the endpoint -- see data-preprocessing
```

### Step 2: Cell Segmentation

```bash
# Mesmer/DeepCell whole-cell (nuclear-first); membrane channels aggregated via the panel column.
# Pass the true acquisition resolution so Mesmer's internal rescaler is correct (~1.0 um IMC).
steinbock segment deepcell --minmax -o masks

# Alternative: Cellpose container (current default model cpsam; channel order reversed vs native)
steinbock segment cellpose --minmax -o masks
```

### Step 3: Single-cell Quantification

```bash
# Extract per-cell MEAN intensities (mean is the default and the right phenotyping aggregator;
# sum confounds cell size with expression)
steinbock measure intensities -o intensities

# Measure cell properties (area, centroid, eccentricity)
steinbock measure regionprops -o regionprops

# Build the spatial neighbor graph (expansion within a max distance; match the graph to the
# biological claim -- contact vs proximity -- in spatial-analysis)
steinbock measure neighbors --type expansion --dmax 15 -o neighbors
```

## Complete Python Workflow

```python
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import squidpy as sq
from pathlib import Path

# === 1. LOAD DATA ===
data_dir = Path('steinbock_output')

intensities = pd.read_csv(data_dir / 'intensities.csv', index_col=0)
regionprops = pd.read_csv(data_dir / 'regionprops.csv', index_col=0)
neighbors = pd.read_csv(data_dir / 'neighbors.csv')

print(f'Loaded {len(intensities)} cells')

# === 2. CREATE ANNDATA ===
adata = ad.AnnData(X=intensities.values, obs=regionprops, var=pd.DataFrame(index=intensities.columns))
adata.obs['image_id'] = [idx.split('_')[0] for idx in intensities.index]
adata.obs['cell_id'] = intensities.index

# Add spatial coordinates
adata.obsm['spatial'] = regionprops[['centroid_y', 'centroid_x']].values

# === 3. PREPROCESSING ===
# Arcsinh transform: cofactor 1 for IMC single-cell means (Hunter 2024), NOT the
# suspension-CyTOF cofactor 5, which over-compresses IMC's lower-count means
adata.layers['counts'] = adata.X.copy()
adata.X = np.arcsinh(adata.X / 1)

# Scale for clustering
sc.pp.scale(adata, max_value=10)
adata.raw = adata.copy()

# === 4. DIMENSIONALITY REDUCTION ===
sc.pp.pca(adata, n_comps=20)
sc.pp.neighbors(adata, n_neighbors=15)
sc.tl.umap(adata)

# === 5. CLUSTERING ===
sc.tl.leiden(adata, resolution=0.8)
print(f'Found {adata.obs["leiden"].nunique()} clusters')

# === 6. PHENOTYPING ===
# Marker expression per cluster
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
marker_genes = sc.get.rank_genes_groups_df(adata, group=None)

# Annotate clusters based on markers
cluster_annotations = {
    '0': 'T cells',
    '1': 'Macrophages',
    '2': 'Tumor',
    '3': 'B cells',
    '4': 'Stromal'
}
adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_annotations)

# === 7. SPATIAL ANALYSIS ===
# Build spatial graph
sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True)

# Neighborhood enrichment
sq.gr.nhood_enrichment(adata, cluster_key='cell_type')

# Co-occurrence analysis
sq.gr.co_occurrence(adata, cluster_key='cell_type')

# Ripley's statistics
sq.gr.ripley(adata, cluster_key='cell_type', mode='L')

# === 8. VISUALIZATION ===
import matplotlib.pyplot as plt

# UMAP by cell type
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
sc.pl.umap(adata, color='cell_type', ax=axes[0], show=False)
sc.pl.umap(adata, color='leiden', ax=axes[1], show=False)
plt.savefig('umap_celltypes.png', dpi=150, bbox_inches='tight')

# Spatial plot
fig, ax = plt.subplots(figsize=(10, 10))
sq.pl.spatial_scatter(adata[adata.obs['image_id'] == 'image1'],
                      color='cell_type', shape=None, size=10, ax=ax)
plt.savefig('spatial_celltypes.png', dpi=150, bbox_inches='tight')

# Neighborhood enrichment heatmap
sq.pl.nhood_enrichment(adata, cluster_key='cell_type')
plt.savefig('neighborhood_enrichment.png', dpi=150, bbox_inches='tight')

# === 9. DIFFERENTIAL ANALYSIS (patient is the unit, NOT the cell) ===
import statsmodels.formula.api as smf

# aggregate to per-image proportions, then test across PATIENTS -- a cell-level or per-image
# test over correlated cells is pseudoreplication and reports p~0 for trivial effects.
# obs must carry patient and condition columns; see differential-analysis for scCODA
# (compositional) and the spatial differential path.
counts = adata.obs.groupby(['patient', 'condition', 'image_id', 'cell_type']).size().unstack(fill_value=0)
image_prop = counts.div(counts.sum(axis=1), axis=0).reset_index()
target = 'Tumor'   # an actual cell_type column from cluster_annotations above (single-word for the formula)
res = smf.mixedlm(f'{target} ~ condition', image_prop, groups=image_prop['patient']).fit()  # patient random effect
print(res.summary())

adata.write('imc_analysis.h5ad')
print('Analysis complete!')
```

## R Alternative (imcRtools)

```r
library(imcRtools)
library(cytomapper)
library(CATALYST)

# Read steinbock output
spe <- read_steinbock('steinbock_output/')

# Transform (cofactor 1 for IMC single-cell means, not 5)
assay(spe, 'exprs') <- asinh(counts(spe) / 1)

# Cluster
spe <- runDR(spe, features = rownames(spe), exprs_values = 'exprs', dr = 'UMAP')
spe <- cluster(spe, features = rownames(spe), exprs_values = 'exprs',
               xdim = 10, ydim = 10, maxK = 20)

# Spatial analysis
spe <- buildSpatialGraph(spe, img_id = 'image_id', type = 'expansion', threshold = 20)
spe <- aggregateNeighbors(spe, colPairName = 'neighborhood', by = 'cluster_id')

# Spatial context
cn <- detectCommunity(spe, colPairName = 'neighborhood',
                       size_threshold = 10, group_by = 'image_id')

# Plot
plotSpatial(spe, img_id = 'image1', node_color_by = 'cluster_id')
```

## QC Checkpoints

| Stage | Check | Action if Failed |
|-------|-------|------------------|
| Preprocessing | No hot pixel streaks | Lower threshold |
| Segmentation | >80% cells detected | Adjust diameter |
| Quantification | All markers extracted | Check panel.csv |
| Clustering | 5-20 clusters | Adjust resolution |
| Spatial | Neighbors detected | Check distance |

## Workflow Variants

### High-plex Panels (40+ markers)
```python
# Use batch-aware clustering
import scvi

scvi.model.SCVI.setup_anndata(adata, batch_key='image_id')
model = scvi.model.SCVI(adata)
model.train()
adata.obsm['X_scvi'] = model.get_latent_representation()
sc.pp.neighbors(adata, use_rep='X_scvi')
```

### Tumor Microenvironment Analysis
```python
# Spatial cell-cell co-location around tumor (per-image, then aggregate to patient).
# Note: sq.gr.ligrec keys ligand-receptor pairs on gene symbols from OmniPath, so it is
# usually empty on a ~40-marker antibody panel -- prefer neighborhood enrichment for IMC.
sq.gr.nhood_enrichment(adata, cluster_key='cell_type')   # see spatial-analysis for the null caveat
```

## Related Skills

- imaging-mass-cytometry/data-preprocessing - Hot pixel, spillover
- imaging-mass-cytometry/cell-segmentation - Cellpose/Mesmer details
- imaging-mass-cytometry/phenotyping - Cluster annotation
- imaging-mass-cytometry/spatial-analysis - Spatial statistics
- imaging-mass-cytometry/differential-analysis - Patient-level cross-condition testing
- imaging-mass-cytometry/interactive-annotation - Manual cell labeling
- imaging-mass-cytometry/quality-metrics - QC metrics
- single-cell/clustering - Clustering methods
- spatial-transcriptomics/spatial-statistics - Related spatial methods
