# Heatmaps and Hierarchical Clustering - Usage Guide

## Overview

Heatmaps render a features-by-samples matrix as a colored grid with optional dendrograms. The choices that determine the biological story are scaling (row z-score vs raw vs robust), distance metric (Euclidean vs correlation vs Manhattan), linkage criterion (ward.D2 vs complete vs average), leaf ordering (default vs optimal), and color mapping (sequential vs diverging, robust bounds). ComplexHeatmap is the modern reference implementation; pheatmap remains maintained for simpler cases; seaborn.clustermap covers Python.

## Prerequisites

```r
install.packages(c('pheatmap', 'seriation', 'dendextend', 'circlize'))
BiocManager::install(c('ComplexHeatmap'))
```

```bash
pip install seaborn scipy pandas numpy
```

## Quick Start

Tell your AI agent what you want to do:
- "Make a clustered heatmap of my expression matrix with row z-score and ward.D2"
- "Add column annotation tracks for condition and batch"
- "Cluster columns within each treatment group, keeping treatment as the major split"
- "Use a robust color scale clipped at the 1st-99th percentile"
- "Render with optimal leaf ordering so visually adjacent rows are actually similar"
- "Save as PDF with raster rendering for a 5000-row matrix"

## Example Prompts

### Bulk RNA-seq expression heatmap

> "Build a publication-quality ComplexHeatmap of the top 500 variable genes, row-z-scored, with ward.D2 linkage on Euclidean distance. Add column annotation for condition and batch. Use a diverging Crameri 'vik' palette with symmetric bounds at the 1st-99th percentile."

### Mutation matrix (OncoPrint)

> "Make an OncoPrint of the recurrent gene mutations across the cohort, color-coded by alteration type (missense, truncating, copy gain, copy loss), sorted by mutation frequency."

### Single-cell pseudobulk

> "Aggregate the scRNA matrix to pseudobulk by cluster, then plot a heatmap of marker genes z-scored across clusters."

### Sample correlation QC

> "Compute pairwise sample correlations on the top-variable-gene matrix; plot as a symmetric heatmap to check for outlier samples and batch grouping."

### Multi-heatmap concatenation

> "Side-by-side heatmap of expression (z-score) and methylation (beta) for the same gene set, sharing the row order from expression clustering."

## What the Agent Will Do

1. Load the input matrix (expression, methylation, peak counts, or generic features-by-samples) and the metadata frame.
2. Filter to the analysis subset - top-variable rows, DE hits, or a pre-specified gene list.
3. Decide scaling: row z-score for cross-sample expression comparison; raw for absolute-interpretable data; robust for outlier-prone data.
4. Decide distance: Euclidean for z-scored data; correlation for co-regulation discovery; binary for mutation matrices.
5. Decide linkage: ward.D2 (default for most genomics); average for WGCNA-style co-expression; complete only when justified.
6. Compute robust color bounds (1st-99th percentile of |matrix|); use symmetric bounds for diverging data.
7. Optionally compute Optimal Leaf Ordering via seriation; pass dendrogram to ComplexHeatmap.
8. Build column and row annotation tracks with explicit color lists.
9. Render with `draw()` (NEVER bare `Heatmap()` in a script); set `use_raster=TRUE, raster_quality=5` for matrices >2000 rows.
10. Export PDF with `cairo_pdf` for vector text + raster cells.

## Tips

- **Always specify `ward.D2` explicitly.** `'ward'` and `'ward.D'` are legacy implementations that do not implement Ward's minimum-variance criterion (Murtagh-Legendre 2014). The names look interchangeable; the dendrograms are not.

- **Use `draw()` in scripts.** A bare `Heatmap(mat)` produces no output and no error inside loops, functions, Quarto chunks, or `Rscript`. Always `draw(Heatmap(mat, ...), merge_legends = TRUE)`.

- **Robust symmetric color bounds.** `bounds <- quantile(abs(mat), 0.99); colorRamp2(c(-bounds, 0, bounds), c('#0072B2', 'white', '#D55E00'))`. Without quantile clipping, one outlier washes out the entire heatmap.

- **Do NOT cluster ordered conditions.** Time points, dose levels, and treatment stages must keep their order. Set `cluster_columns = FALSE`; use `column_split` to group while preserving order.

- **`z_score=0` vs `standard_scale=0` in seaborn are different.** `z_score=0` standardizes rows to mean 0 / sd 1. `standard_scale=0` rescales rows to [0, 1] - looks similar but compresses outliers non-linearly. The two are mutually exclusive.

- **`gaps_col` is ignored when `cluster_cols=TRUE`** in pheatmap. Either disable clustering or switch to ComplexHeatmap `column_split`.

- **Optimal Leaf Ordering (OLO)** rearranges leaves within a fixed dendrogram structure to minimize adjacent-leaf dissimilarity (Bar-Joseph 2001). Reveals block structure the default ordering hides. Use `seriation::seriate(d, method='OLO')`.

- **Correlation distance preserves co-regulation but is sensitive to batch.** Always batch-correct before clustering with correlation distance, or modules cluster by batch.

- **Sparse matrices break z-scoring.** Rows with most-zero values have unreliable sd. Filter low-expression rows before z-scoring, or use robust scaling (median/mad).

- **`raster_quality = 5` for publication.** The default raster quality looks pixelated when scaled for print.

- **pheatmap is still maintained.** v1.0.13 as of 2025-06. Not abandoned. For simple cases it is sufficient; for split panels, multi-heatmap concatenation, or oncoPrint use ComplexHeatmap.

- **Single-cell pseudobulk before heatmap.** Plotting per-cell expression for >10000 cells crashes PDF renderers. Aggregate to cluster × gene pseudobulk first.

## Related Skills

- data-visualization/color-palettes - Sequential and diverging palette selection
- data-visualization/oncoprint-mutation-matrices - Specialized mutation-matrix heatmap
- data-visualization/multipanel-figures - Combining heatmaps into journal figures
- single-cell/markers-annotation - scRNA dotplot / matrixplot alternatives
- differential-expression/de-visualization - DE-gene heatmap recipes
