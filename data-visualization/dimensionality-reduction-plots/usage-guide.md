# Dimensionality-Reduction Plots - Usage Guide

## Overview

PCA, t-SNE, UMAP, and PHATE are the four standard methods for projecting high-dimensional omics data to 2D. Each preserves a different property - variance (PCA), local neighborhoods (t-SNE), local + partial global (UMAP), continuous transitions (PHATE). The Chari-Pachter 2023 critique established that 2D embeddings lose >95% of high-dimensional geometry, so embeddings communicate "these points are similar locally" and nothing more. Hyperparameters, random seeds, and explicit interpretation limits matter for reproducibility.

## Prerequisites

```bash
pip install scanpy umap-learn openTSNE phate scikit-learn matplotlib
```

```r
install.packages(c('Rtsne', 'uwot', 'PCAtools', 'phateR'))
BiocManager::install(c('PCAtools', 'DESeq2'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Make a PCA plot of bulk RNA-seq for sample QC, colored by condition and shaped by batch"
- "Compute UMAP from scanpy AnnData using n_neighbors=30, min_dist=0.3, random_state=42"
- "Run t-SNE with Kobak-Berens defaults: PCA init, perplexity 30, learning_rate n/12"
- "Plot PHATE for trajectory display instead of UMAP"
- "Annotate PCA axes with variance explained percentages"
- "Show loadings as arrows on PC1 vs PC2"

## Example Prompts

### Bulk PCA for sample QC

> "PCA on vst-normalized bulk RNA-seq counts. Plot PC1 vs PC2 colored by condition, shaped by batch. Axis labels must include variance explained. Add screeplot of PC1-10."

### Single-cell UMAP

> "From scanpy AnnData, compute PCA(50), neighbors(30, n_pcs=50), UMAP(min_dist=0.3, random_state=42). Plot colored by Leiden cluster with on-data labels."

### Kobak-Berens t-SNE

> "Run openTSNE with PCA initialization, perplexity=30, learning_rate=n/12. Plot the embedding colored by cluster assignment."

### PHATE for trajectory

> "Use PHATE instead of UMAP for displaying developmental data - preserves continuous transitions."

### Method comparison

> "Run PCA, UMAP, t-SNE, and PHATE on the same matrix. Display side-by-side. Annotate each panel with its preservation property."

## What the Agent Will Do

1. Decide method: PCA for variance-explained QC; UMAP for cluster overview; t-SNE if cluster boundaries are critical; PHATE for continuous trajectories.
2. Pre-process input: normalize, log-transform, scale (for PCA); compute PCA(50) before t-SNE/UMAP for single-cell.
3. Set explicit hyperparameters: perplexity / n_neighbors / min_dist / random_state.
4. Fit the projection with a fixed seed for reproducibility.
5. Plot with axis labels: PCA shows variance %; UMAP/t-SNE/PHATE label only "UMAP1 / UMAP2" (no units).
6. Color by categorical (CVD-safe palette) or continuous (perceptually-uniform colormap).
7. Annotate with cluster labels on-plot OR via legend depending on cluster count.
8. State in caption: hyperparameters used; embedding's interpretation limit.

## Tips

- **2D embeddings lose >95% of high-dim geometry** (Chari-Pachter 2023). Distance between clusters is NOT meaningful. Density within clusters is dominated by `min_dist`, not biology. State this limit in the caption.

- **Always set `random_state`** (umap-learn, sklearn) or `seed=` (uwot, Rtsne). Without it, layouts vary across runs.

- **For t-SNE, use Kobak-Berens defaults**: `init='pca'`, perplexity 30, learning_rate = n/12. Default learning_rate=200 over-shrinks large data.

- **Always label PCA axes with variance explained**: `PC1 (45%) PC2 (12%)`. Without this, "clusters" at PC1=4%, PC2=3% may be noise.

- **PC1 = library size is the canonical pre-processing failure.** Use `vst()` / `rlog()` (DESeq2) or `log + scale` before PCA.

- **`n_neighbors` controls local-vs-global trade-off in UMAP.** Small (5-10) = fragmented; large (50+) = merged. 15-30 is standard.

- **`min_dist` controls tightness, NOT cluster membership.** Smaller = denser; does not change which cells cluster together.

- **Do not interpret cluster shape.** UMAP/t-SNE cluster shape is an embedding artifact, not biology. Cluster membership is the biological observation.

- **Loadings only exist for PCA.** UMAP/t-SNE have no loadings; "what gene drives the axis" requires PCA.

- **scanpy save trap**: `sc.pl.umap(save='_x.pdf')` writes to `sc.settings.figdir + 'umap_x.pdf'`, not cwd. Default DPI 150 - set 300+ for publication.

- **Trajectory claims require validation.** Don't infer trajectories from UMAP "gaps." Use RNA velocity, diffusion pseudotime, or PHATE - and validate against orthogonal markers.

- **Batch correction MUST precede UMAP.** UMAP preserves local neighborhoods and can hide batch effects that PCA exposes. Diagnose batch in PCA; correct upstream; then UMAP.

## Related Skills

- single-cell/preprocessing - PCA / neighbor graph before embedding
- single-cell/clustering - Leiden / Louvain assignments to color UMAP
- single-cell/trajectory-inference - Pseudotime / RNA velocity for trajectory claims
- data-visualization/color-palettes - Categorical and perceptual palettes
- data-visualization/distribution-plots - Per-cluster gene-expression follow-up
