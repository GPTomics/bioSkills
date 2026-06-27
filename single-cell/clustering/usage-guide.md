# Single-Cell Clustering - Usage Guide

## Overview

This skill covers dimensionality reduction and graph-based clustering for single-cell RNA-seq using Scanpy (Python) and Seurat (R). It treats clustering as exploratory hypothesis generation: resolution has no biological meaning, clusters must be validated by markers and stability before being named cell populations, and UMAP/tSNE are for display only. The skill resolves the hard decisions - Leiden vs Louvain, how many PCs, how to sweep and validate resolution, when a split is over-clustering, and why post-clustering marker p-values are not valid inference.

## Prerequisites

**Python (Scanpy):**
```bash
pip install scanpy igraph leidenalg matplotlib
```

**R (Seurat):**
```r
install.packages(c('Seurat', 'clustree'))
```

## Quick Start

Tell the agent what to do:
- "Run PCA and cluster my single-cell data"
- "Sweep clustering resolutions and tell me which to use"
- "Are these two clusters really separate cell types?"
- "Generate a UMAP colored by cluster"

## Example Prompts

### Dimensionality Reduction
> "Run PCA and show the elbow plot to pick the number of PCs"

> "How many PCs should drive the neighbor graph for this dataset?"

### Clustering
> "Cluster with Leiden and pin the backend so it is reproducible"

> "Sweep resolutions 0.2 to 1.0 and build a clustree to choose granularity"

> "Subcluster cluster 3 and recompute HVGs and PCA on the subset"

### Validation
> "Check whether my clusters are stable to bootstrapping"

> "Run a significance test to decide if this split is over-clustering"

> "These two clusters share all their top markers - should they be merged?"

### Visualization
> "Show clusters on a UMAP and overlay CD3D, MS4A1, and CD14"

## What the Agent Will Do

1. Run PCA and choose n_pcs from the elbow (n_pcs dominates the result more than n_neighbors).
2. Build the kNN graph and partition with Leiden, pinning the backend (`flavor='igraph', n_iterations=2, directed=False`) for reproducibility.
3. Sweep a resolution range rather than fixing one value, and visualize cell flow with clustree.
4. Validate clusters: stability by bootstrapping, distinct markers per cluster, and a significance test (scSHC/CHOIR) before claiming populations.
5. Adjudicate biological vs technical: rule out cell-cycle, dissociation stress, mitochondrial, ambient, and batch drivers.
6. Generate UMAP/tSNE for display, treating distances as non-metric.

## Tips

- **Leiden over Louvain** - Louvain can return disconnected communities (Traag 2019); Scanpy uses Leiden, Seurat defaults to Louvain (`algorithm=1`).
- **Pin the Leiden backend** - Scanpy 1.10 is mid-migration; unpinned `flavor`/`n_iterations` gives different labels across versions.
- **n_pcs is the big lever** - it matters far more than n_neighbors; set it from the elbow, then tune neighbors.
- **Resolution is not a truth knob** - sweep it, use clustree, pick the coarsest defensible level; tuning to match a reference is confirmation bias.
- **Stable does not mean real** - a reproducible cluster can be cell-cycle, stress, mito, ambient, or batch; pair stability with a significance test.
- **Over-clustering is the default failure** - if adjacent clusters share all markers, merge them.
- **Marker p-values are not inference** - clustering then testing on the same data is double-dipping; use markers for ranking, validate with scSHC/CHOIR + ClusterDE.
- **UMAP is for the eye, not the ruler** - cluster on the graph; never read embedding distances as biology (Chari & Pachter 2023).

## Related Skills

- preprocessing - QC, normalization, and HVG selection that must precede clustering
- doublet-detection - Remove doublets before clustering so they do not form fake intermediate clusters
- batch-integration - Integrate batches before clustering when a cluster tracks a single sample
- markers-annotation - Find and test markers per cluster (with the double-dipping caveat)
- cell-annotation - Assign cell-type identities to validated clusters
- single-cell/differential-abundance - Test whether cluster proportions shift across conditions
- data-visualization/dimensionality-reduction-plots - Publication-quality UMAP/tSNE/PCA figures
- pathway-analysis/go-enrichment - Interpret per-cluster marker sets
