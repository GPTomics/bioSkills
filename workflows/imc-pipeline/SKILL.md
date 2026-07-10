---
name: bio-workflows-imc-pipeline
description: Orchestrates imaging mass cytometry from raw MCD acquisitions to patient-level spatial analysis, chaining steinbock preprocessing, Mesmer/Cellpose segmentation, single-cell quantification, phenotyping, and squidpy spatial statistics. Use when committing the panel + segmentation frame + pixel size (every per-cell number is a mask-bounded pixel average), compensating channel spillover on PIXELS before segmentation but running REDSEA lateral-spillover on the per-cell table AFTER segmentation, using arcsinh cofactor 1 (not the suspension-CyTOF 5), and aggregating to the PATIENT before any cross-condition test (cells and ROIs from one patient are not independent replicates). Hands mechanism to the imaging-mass-cytometry component skills; not a re-teach of any single step.
tool_type: python
primary_tool: steinbock
goal_approach_exempt: true
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

Reference examples tested with: Cellpose 4.0+ (cpsam model), anndata 0.10+, matplotlib 3.8+, numpy 1.26+, pandas 2.2+, scanpy 1.10+, scvi-tools 1.1+, squidpy 1.3+, steinbock 0.16+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Imaging Mass Cytometry Pipeline

**"Process my imaging mass cytometry data from images to spatial analysis"** -> Orchestrate image preprocessing (steinbock), cell segmentation (Cellpose), phenotyping (FlowSOM/scanpy), spatial neighborhood analysis (squidpy), and tissue community detection.

This is a workflow skill: it owns the chaining decisions and hand-offs, not the internals of any one step.

## The governing principle

Segmentation is the largest irreversible error source, and it is spatial: every per-cell number is a mask-bounded pixel average, so a wrong boundary fabricates cell types before any expression QC can see them. The seam ORDER — and the patient-level unit — is therefore what decides trustworthiness.

1. **The comparison frame (panel + segmentation frame + pixel size) is committed once and inherited by every per-cell number.** The summed membrane channel encodes a cell-type bias (sum BROADLY-expressed markers, or segmentation under-performs on cell types lacking a strong membrane marker); Mesmer was trained at model_mpp ~0.5 um and rescales the input to it, so passing the wrong pixel size (Mesmer's `image_mpp` defaults to None = NO rescaling, assuming the input is already at model resolution — the true pixel size must be passed explicitly; steinbock's `--pixelsize` flag wraps `image_mpp` and defaults to 1.0) rescales cells to the wrong learned size and degrades every boundary. No downstream step recovers a merged or split cell.
2. **Channel spillover is compensated on PIXELS before segmentation; lateral spillover (REDSEA) runs on the per-cell table AFTER segmentation — they are DIFFERENT problems.** Metal-isotope crosstalk is a pixel-level NNLS correction whose compensated value must be what gets averaged into the per-cell mean (post-aggregation is wrong). REDSEA corrects real signal leaking across shared cell boundaries at ~1 um even with perfect segmentation and zero channel spillover — it is defined on segmented neighbors, so it must run after segmentation. Running REDSEA pre-segmentation, or channel comp post-aggregation, is a category error.
3. **The experimental unit is the PATIENT, not the cell or the ROI.** Cells and ROIs from one patient are not independent replicates; a cell-level or per-image test over correlated cells is pseudoreplication (reports p~0 for trivial effects). Aggregate to per-patient proportions/summaries, then a mixed model / scCODA. Arcsinh cofactor is 1 for IMC integer ion counts, NOT the suspension-CyTOF 5 (which over-compresses them). Impossible lineage-exclusive co-expression is a segmentation/spillover ALARM, not a hybrid cell type.

## Made-once commitments

| Commitment | Consequence inherited downstream |
|------------|----------------------------------|
| Panel (metal->antibody; membrane-sum channels) | Which channels extract and phenotype; a narrow membrane sum biases segmentation against some cell types |
| Segmentation frame (nuclear + membrane channels) | Every per-cell number (all are mask-bounded pixel averages); the largest irreversible error source |
| Pixel size (steinbock `--pixelsize` / Mesmer `image_mpp`, ~1.0 um for IMC) | Boundary quality + all spatial distances; the wrong value rescales cells to the wrong learned size |
| Arcsinh cofactor = 1 (IMC), not 5 (CyTOF) | Clustering/phenotyping distances; cofactor 5 over-compresses integer ion counts |

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
# --pixelsize is steinbock's CLI flag for the acquisition resolution (it wraps Mesmer's image_mpp);
# steinbock defaults it to 1.0 um for IMC, so pass the true value explicitly rather than relying on it.
steinbock segment deepcell --pixelsize 1.0 --minmax -o masks

# Alternative: Cellpose container (Cellpose 4+ default model cpsam; channel order reversed vs native)
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
adata.obs['image_id'] = pd.Categorical([idx.rsplit('_', 1)[0] for idx in intensities.index])   # strip only the trailing cell index: rsplit keeps Patient1_ROI002 distinct from Patient1_ROI001. squidpy library_key requires a categorical, not object/string
adata.obs['cell_id'] = intensities.index

# Add spatial coordinates (skimage regionprops_table names them centroid-0 (y) / centroid-1 (x))
adata.obsm['spatial'] = regionprops[['centroid-0', 'centroid-1']].values

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
# Build spatial graph PER IMAGE (library_key), else Delaunay fabricates edges across ROIs
sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True, library_key='image_id')

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

# Spatial plot. Pick the image dynamically: image_id is derived from the cell index, so a hardcoded
# literal selects zero cells and spatial_scatter errors on the empty subset.
fig, ax = plt.subplots(figsize=(10, 10))
first_image = adata.obs['image_id'].iloc[0]
sq.pl.spatial_scatter(adata[adata.obs['image_id'] == first_image],
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
counts = adata.obs.groupby(['patient', 'condition', 'image_id', 'cell_type'], observed=True).size().unstack(fill_value=0)   # observed=True: image_id is categorical; the default expands the full cartesian product into all-zero phantom rows -> NaN proportions
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

# Cluster (CATALYST runDR takes assay=; cluster() always uses the 'exprs' assay, no assay arg)
spe <- runDR(spe, features = rownames(spe), assay = 'exprs', dr = 'UMAP')
spe <- cluster(spe, features = rownames(spe), xdim = 10, ydim = 10, maxK = 20)

# Spatial analysis. buildSpatialGraph names the colPair '<type>_interaction_graph';
# aggregateNeighbors counts a label via aggregate_by='metadata' + count_by=.
spe <- buildSpatialGraph(spe, img_id = 'sample_id', type = 'expansion', threshold = 20)
spe <- aggregateNeighbors(spe, colPairName = 'expansion_interaction_graph',
                          aggregate_by = 'metadata', count_by = 'cluster_id')

# Spatial context
spe <- detectCommunity(spe, colPairName = 'expansion_interaction_graph',
                       size_threshold = 10, group_by = 'sample_id')

# Plot (img_id is the colData COLUMN used to facet; read_steinbock names it 'sample_id', not 'image_id')
plotSpatial(spe, img_id = 'sample_id', node_color_by = 'cluster_id')
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

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Impossible double-positive "hybrid" cell types | Spillover not corrected before phenotyping (channel and/or lateral) | NNLS channel compensation on pixels before segmentation; REDSEA on the per-cell table after; treat lineage-exclusive co-expression as a QC failure until proven |
| Every boundary degraded, cells the wrong size | Wrong pixel size (Mesmer `image_mpp` defaults None=no rescaling, model trained at ~0.5; steinbock `--pixelsize` defaults 1.0) | Pass the true acquisition resolution explicitly (~1.0 um for IMC) |
| Macrophages under-captured; biased comparison | Nuclear-expansion segmentation cross-compared with whole-cell data | Never quantitatively compare expansion-segmented vs whole-cell; report the expansion radius; use constrained (not free) dilation |
| p~0 for a trivial effect | Pseudoreplication (cells/ROIs treated as replicates) | Aggregate to per-patient summaries; mixed model with patient random effect / scCODA |
| Markers over-compressed, noise clusters | Arcsinh cofactor 5 used on IMC | Cofactor 1 for IMC integer ion counts |
| Acquisition batch drives the clusters | Batch confounded with / not modeled against condition | Randomize acquisition order; batch-aware clustering (Harmony/scVI) for clustering ONLY; model batch as a covariate; no rescue if batch==condition |

## References

- Windhager J, Zanotelli VRT, Schulz D, et al (2023) An end-to-end workflow for multiplexed image processing and analysis. *Nature Protocols* 18:3565-3613. DOI 10.1038/s41596-023-00881-0. (steinbock.)
- Greenwald NF, Miller G, Moen E, et al (2022) Whole-cell segmentation of tissue images with human-level performance using large-scale data annotation and deep learning. *Nature Biotechnology* 40:555-565. DOI 10.1038/s41587-021-01094-0. (Mesmer/DeepCell.)
- Bai Y, Zhu B, Rovira-Clave X, et al (2021) Adjacent cell marker lateral spillover compensation and reinforcement for multiplexed images. *Frontiers in Immunology* 12:652631. DOI 10.3389/fimmu.2021.652631. (REDSEA.)
- Palla G, Spitzer H, Klein M, et al (2022) Squidpy: a scalable framework for spatial omics analysis. *Nature Methods* 19:171-178. DOI 10.1038/s41592-021-01358-2.
- Hunter B, Nicorescu I, Foster E, et al (2024) OPTIMAL: an OPTimized Imaging Mass cytometry AnaLysis framework for benchmarking segmentation and data exploration. *Cytometry Part A* 105:36-53. DOI 10.1002/cyto.a.24803. (arcsinh cofactor 1 for IMC.)

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
