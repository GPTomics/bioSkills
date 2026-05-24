---
name: bio-data-visualization-heatmaps-clustering
description: Build clustered heatmaps for expression matrices and other features-by-samples data with rigorous distance/linkage/scaling choices, robust color mapping, optimal leaf ordering, and ComplexHeatmap/pheatmap/seaborn rendering. Covers the ward.D vs ward.D2 trap, the row-vs-column scaling decision, multi-track annotations, oncoPrint, and raster rendering for large matrices. Use when visualizing expression patterns across samples or identifying co-regulated clusters.
tool_type: mixed
primary_tool: ComplexHeatmap
---

## Version Compatibility

Reference examples tested with: ComplexHeatmap 2.18+, pheatmap 1.0.13 (still maintained as of 2025-06), circlize 0.4.16+, seaborn 0.13+, scipy 1.12+, scanpy 1.10+, ggplot2 3.5+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Heatmaps and Hierarchical Clustering

**"Make a clustered heatmap"** -> Render an expression / feature matrix as a colored grid with hierarchical-clustering dendrograms, after committing to (a) how to scale the data (row z-score vs raw vs robust), (b) which distance metric (Euclidean vs correlation vs Manhattan), (c) which linkage criterion (ward.D2 vs complete vs average), (d) how to order the leaves (default vs optimal leaf ordering), (e) how to map values to color (sequential vs diverging, robust quantile bounds), and (f) which package can handle the matrix size and annotation complexity.

- R: `ComplexHeatmap::Heatmap()` (modern default), `pheatmap::pheatmap()` (still maintained, simpler API)
- Python: `seaborn.clustermap()`, `scanpy.pl.heatmap()` (single-cell-aware)

## The Single Most Important Modern Insight -- Distance, Linkage, and Scaling Are Three Independent Decisions

A heatmap's dendrograms are produced by three orthogonal choices, each with material biological consequences:

1. **Scaling** decides what "similar" means. Row z-scoring asks "do these genes covary across samples?" — it strips absolute level. Raw values ask "do these genes have similar magnitude AND pattern?" Robust scaling (quantile clip) asks "do these covary after suppressing outliers?"
2. **Distance metric** decides how dissimilar two profiles are. Euclidean on z-scored data ≈ `1 − Pearson` correlation; Manhattan tolerates outliers; correlation distance preserves co-regulation patterns regardless of amplitude.
3. **Linkage** decides how to merge clusters. Ward minimizes within-cluster sum of squares (compact, spherical clusters); complete uses max distance (compact, outlier-sensitive); average is balanced; single chains (almost never what genomics wants).

The biological story changes depending on these choices. A "module" identified with `complete` linkage on Euclidean distance of raw counts is *not the same module* identified with `ward.D2` on correlation distance of z-scored data. Both can be defensible; neither is automatic. **Pick deliberately, document, and verify the clustering against orthogonal evidence before claiming the modules are biological.**

## The ward.D vs ward.D2 Trap (Murtagh-Legendre 2014)

R `stats::hclust` exposes two methods both labeled "Ward": `ward.D` and `ward.D2`. They produce *different* dendrograms on the same data. Only `ward.D2` (squared distances input) implements Ward's actual minimum-variance criterion (Murtagh & Legendre 2014 *J Classif* 31:274). `ward.D` is a historical implementation that does not.

```r
hclust(dist(x), method = 'ward.D')   # NOT Ward's criterion -- legacy
hclust(dist(x), method = 'ward.D2')  # Ward's actual minimum-variance criterion
```

`pheatmap::pheatmap(clustering_method='ward.D2')` and `ComplexHeatmap::Heatmap(clustering_method_rows='ward.D2')` both pass through to `hclust`. Always specify `ward.D2` unless reproducing a paper that used the unlabeled `ward` (which actually called `ward.D` pre-R 3.1).

## Decision Tree by Scenario

| Scenario | Scaling | Distance | Linkage | Why |
|----------|---------|----------|---------|-----|
| Bulk RNA-seq, expression patterns across samples | row z-score | euclidean | ward.D2 | Standard; z-score removes absolute level so co-regulated genes cluster regardless of magnitude |
| Methylation beta values (already bounded [0,1]) | raw (no scale) | euclidean or manhattan | ward.D2 | Beta values are interpretable on absolute scale; scaling would distort |
| Co-expression module discovery | row z-score | correlation (`1 - cor`) | average | WGCNA convention; preserves co-regulation pattern |
| ChIP/ATAC peak intensity across samples | raw log-counts | euclidean | ward.D2 | Peaks are interpretable on absolute scale after log |
| Sample QC (correlation of samples) | column-wise raw | correlation | ward.D2 | The correlation IS the data; don't scale before computing it |
| Methylation array with outliers | raw + clip 1-99% | euclidean | ward.D2 | Outliers dominate Euclidean; robust clip preserves signal |
| Single-cell pseudobulk by cell type | row z-score | euclidean | ward.D2 | Same as bulk; downsample to <500 cells per type for rendering |
| Mutation matrix (binary present/absent) | raw | binary (jaccard) | average or complete | Standard distance for binary data; ward inappropriate |
| Drug response across cell lines | row z-score | spearman correlation | ward.D2 | Drug-rank patterns matter more than absolute IC50 |

## Color Mapping -- The Quietly Most-Important Choice

A heatmap is a color encoding of a matrix. The default linear mapping from data to color is rarely correct:

1. **Diverging data needs symmetric bounds.** For z-scores or log-fold changes, the color bar must be symmetric around zero. `colorRamp2(c(-2, 0, 2), c('#0072B2', 'white', '#D55E00'))` ALWAYS, not `colorRamp2(c(min, mean, max), ...)`.

2. **Robust quantile bounds.** Single outliers compress the entire color scale. Clip at 1st/99th percentile before mapping: `bounds <- quantile(mat, c(0.01, 0.99))`. ComplexHeatmap's `colorRamp2(c(bounds[1], 0, bounds[2]), ...)` is standard. Without this, one outlier sample turns the entire heatmap pale.

3. **Sequential data uses a perceptually-uniform colormap.** viridis, magma, cividis (Nuñez 2018), or batlow (Crameri 2020). NOT jet, NOT rainbow, NOT `colorRampPalette(c('blue','red'))(100)` which has a non-monotonic luminance.

4. **Diverging palettes from Crameri** (`vik`, `roma`) or ColorBrewer (`RdBu`, `BrBG`) are perceptually uniform. Reverse the default direction for log-fold-change (negative = blue, positive = red, by biological convention).

## Optimal Leaf Ordering (Bar-Joseph 2001)

A dendrogram for n leaves has 2^(n-1) consistent linear orderings — only one is the leaf order shown. Default `hclust` gives a deterministic but visually arbitrary ordering. **Optimal Leaf Ordering (OLO)** chooses the consistent ordering that minimizes the sum of distances between adjacent leaves — making visually adjacent rows actually similar, and revealing block structure in the heatmap that the default ordering hides.

```r
library(ComplexHeatmap)
library(seriation)

# OLO via seriation
dist_rows <- dist(mat)
hc_rows <- hclust(dist_rows, method = 'ward.D2')
olo_rows <- seriate(dist_rows, method = 'OLO', control = list(hclust = hc_rows))

Heatmap(mat,
        cluster_rows = as.dendrogram(olo_rows[[1]]),
        cluster_columns = TRUE,
        clustering_method_columns = 'ward.D2')
```

For matrices >2000 rows OLO becomes slow (O(n^4) in worst case; modern implementations are much faster). The trade-off is worth it for publication figures.

## Annotation Tracks -- ComplexHeatmap as the Reference

**Goal:** Render an annotated heatmap with column metadata (condition, batch, age), row metadata (pathway, gene class), and split panels for grouped display.

**Approach:** Define `HeatmapAnnotation` (column) and `rowAnnotation` objects with explicit color lists; render with `Heatmap()` specifying `row_split`/`column_split` for grouped layout; use `draw()` to commit, not bare `Heatmap()`, when running non-interactively.

```r
library(ComplexHeatmap)
library(circlize)

# Robust symmetric color mapping
bounds <- quantile(abs(mat[!is.na(mat)]), 0.99)
col_fun <- colorRamp2(c(-bounds, 0, bounds), c('#0072B2', 'white', '#D55E00'))

# Column metadata
ha_col <- HeatmapAnnotation(
    Condition = metadata$condition,
    Batch     = metadata$batch,
    Age       = anno_barplot(metadata$age),
    col = list(
        Condition = c(Control = '#56B4E9', Treatment = '#D55E00'),
        Batch     = c(A = '#009E73', B = '#0072B2', C = '#CC79A7')
    ),
    annotation_name_gp = gpar(fontsize = 8)
)

# Row metadata
ha_row <- rowAnnotation(
    Pathway = gene_info$pathway,
    LogFC   = anno_barplot(gene_info$log2FC, baseline = 0,
                            gp = gpar(fill = ifelse(gene_info$log2FC > 0,
                                                     '#D55E00', '#0072B2'))),
    col = list(Pathway = c(Metabolism = '#8491B4', Signaling = '#91D1C2'))
)

ht <- Heatmap(mat,
              name = 'Z-score',
              col  = col_fun,
              top_annotation  = ha_col,
              left_annotation = ha_row,
              row_split    = gene_info$pathway,
              column_split = metadata$condition,
              clustering_method_rows    = 'ward.D2',
              clustering_method_columns = 'ward.D2',
              clustering_distance_rows    = 'euclidean',
              clustering_distance_columns = 'euclidean',
              show_row_names = FALSE,
              use_raster = TRUE)          # rasterize cell layer for >2000 rows

draw(ht, merge_legends = TRUE)            # draw() not bare Heatmap()
```

### The `draw()` requirement (silent failure)

A bare `Heatmap(mat)` works at the R console because auto-print invokes `draw()`. **Inside `for`, `lapply`, `function`, Quarto/Rmd chunks, or `Rscript`, a bare `Heatmap()` produces no output and no error.** Always wrap in `draw()` non-interactively. Only `draw()` exposes `merge_legends`, `heatmap_legend_side`, `ht_gap`, and `padding`.

## seaborn.clustermap (Python)

```python
import seaborn as sns
import numpy as np
import pandas as pd

# Robust symmetric bounds (1-99% quantile)
vmax = np.quantile(np.abs(df.values[~np.isnan(df.values)]), 0.99)

# col_colors / row_colors for categorical annotations
condition_colors = metadata['condition'].map({'Control': '#56B4E9', 'Treatment': '#D55E00'})
batch_colors = metadata['batch'].map({'A': '#009E73', 'B': '#0072B2', 'C': '#CC79A7'})
col_colors = pd.DataFrame({'Condition': condition_colors, 'Batch': batch_colors})

g = sns.clustermap(df,
                   cmap='RdBu_r', center=0, vmin=-vmax, vmax=vmax,
                   row_cluster=True, col_cluster=True,
                   method='ward',                      # seaborn uses scipy ward, equivalent to R ward.D2
                   metric='euclidean',
                   z_score=0,                          # 0 = rows, 1 = columns
                   col_colors=col_colors,
                   dendrogram_ratio=0.15,
                   cbar_pos=(0.02, 0.8, 0.03, 0.15),
                   figsize=(10, 12),
                   rasterized=True)                    # rasterize the cell layer
```

**`standard_scale` vs `z_score` confusion:**
- `z_score=0` standardizes ROWS to mean 0, SD 1 (most common request)
- `z_score=1` standardizes COLUMNS to mean 0, SD 1
- `standard_scale=0` rescales ROWS to [0, 1] via `(x − min) / (max − min)` — NOT z-scoring, compresses outliers nonlinearly
- The two are mutually exclusive — passing both errors

A heatmap published with `standard_scale` looks like a z-scored heatmap but the color encoding is not interpretable as standard deviations.

## OncoPrint -- The Specialized Mutation-Matrix Heatmap

OncoPrint (Cerami 2012 *Cancer Discov* 2:401; canonical at cBioPortal) is a stylized heatmap for mutation matrices where each cell encodes multiple alteration types via overlapping rectangles. Different from generic heatmaps — see `data-visualization/oncoprint-mutation-matrices` for the dedicated skill. Mentioned here only to note that `ComplexHeatmap::oncoPrint()` is the R implementation and inherits all the cluster/annotation machinery of `Heatmap()`.

## Per-Method Failure Modes

### ward.D used when ward.D2 was intended

**Trigger:** `clustering_method = 'ward'` or `'ward.D'` (with or without the .D).

**Mechanism:** R `hclust` `ward.D` is a legacy implementation that does NOT use squared distances — it does not implement Ward's minimum-variance criterion (Murtagh-Legendre 2014).

**Symptom:** Different dendrogram than published papers that say "Ward"; clusters look subtly different; reproducibility issues across R versions.

**Fix:** Always specify `ward.D2` explicitly. For pheatmap and ComplexHeatmap, pass `clustering_method_rows = 'ward.D2'` and `clustering_method_columns = 'ward.D2'`.

### One outlier compresses the color scale

**Trigger:** Plotting matrix without quantile clipping; one extreme value dominates the color range.

**Mechanism:** Default `colorRamp2(c(min(mat), 0, max(mat)), ...)` is dominated by the outlier — the rest of the matrix renders within a narrow band of pale colors.

**Symptom:** Heatmap looks "washed out" except for one cell or one column; biological pattern invisible.

**Fix:** `bounds <- quantile(abs(mat), 0.99); col_fun <- colorRamp2(c(-bounds, 0, bounds), c('#0072B2', 'white', '#D55E00'))`. ComplexHeatmap's `colorRamp2` does NOT clip values exceeding the range — they render at the extreme color, which is the intended behavior.

### ComplexHeatmap silently produces no output in a script

**Trigger:** Bare `Heatmap(mat)` inside a `for` loop, `lapply`, `function()`, Quarto/Rmd code chunk, or `Rscript` invocation.

**Mechanism:** Auto-print only happens at the top-level R prompt; in non-interactive contexts the Heatmap object is created but never rendered.

**Symptom:** No error, no warning, no PDF output. Looks like the script ran successfully.

**Fix:** Always wrap in `draw()`: `pdf('out.pdf', ...); draw(Heatmap(mat, ...)); dev.off()`. Use the `draw()` call to set legend layout: `draw(ht, merge_legends = TRUE, heatmap_legend_side = 'right')`.

### Clustering applied to ordered conditions

**Trigger:** `cluster_columns = TRUE` when columns are an ordered sequence (time points, dose levels, treatment stages).

**Mechanism:** Hierarchical clustering re-orders columns to maximize within-cluster similarity, destroying the time/dose axis.

**Symptom:** Time-course heatmap with time points scrambled; reader cannot follow temporal pattern.

**Fix:** `cluster_columns = FALSE` for ordered conditions. To group while preserving order, use `column_split` or `column_order` explicitly.

### Z-score on a sparse matrix

**Trigger:** Row z-scoring a matrix with many zero values (e.g., single-cell expression, sparse peak counts).

**Mechanism:** `(x − mean) / sd` is ill-defined when a row is mostly zeros — sd is dominated by the few non-zero values; z-scores explode for the non-zero entries.

**Symptom:** A few cells render as extreme colors; most cells are washed-out near-zero.

**Fix:** For single-cell, work with cluster-summarized pseudobulk matrices, not raw single-cell expression. For sparse peak data, filter rows by minimum non-zero count before scaling.

### Correlation distance applied to data with batch effect

**Trigger:** `clustering_distance_rows = 'correlation'` on a matrix where samples have a strong batch effect.

**Mechanism:** Correlation preserves co-regulation pattern but is sensitive to global structure. If batch shifts ALL genes up in one batch, correlation distance reads this as "co-regulation."

**Symptom:** Modules cluster by batch, not by biology.

**Fix:** Batch-correct (limma::removeBatchEffect or ComBat) before clustering. Confirm with PCA that batch is no longer the dominant axis.

### pheatmap's `gaps_col` interacts with `cluster_cols`

**Trigger:** Setting `gaps_col = c(5, 10)` with `cluster_cols = TRUE`.

**Mechanism:** When `cluster_cols = TRUE`, `gaps_col` is silently ignored; the dendrogram-determined order has no concept of position.

**Symptom:** No gaps appear; no warning.

**Fix:** To use `gaps_col`, set `cluster_cols = FALSE` and pre-arrange the columns explicitly. Or in ComplexHeatmap use `column_split` to combine clustering AND visual gaps.

### Raster rendering at low DPI looks pixelated

**Trigger:** `use_raster = TRUE` (default for large heatmaps in ComplexHeatmap) with default `raster_quality = 1`.

**Mechanism:** Default raster quality is set for screen rendering; the bitmap is upscaled for PDF output, producing blocky cells.

**Symptom:** Heatmap cells look pixelated at print zoom; diagonal "stair-stepping" on cell boundaries.

**Fix:** `Heatmap(..., use_raster = TRUE, raster_quality = 5)` increases the raster resolution. Set `raster_device = 'CairoPNG'` for transparency support.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| ComplexHeatmap and pheatmap give different dendrograms | pheatmap uses `dist()` with default `method = 'euclidean'`; ComplexHeatmap defaults to the same but different clustering distance defaults | Verify both with `clustering_distance_rows = 'euclidean'`, `clustering_method_rows = 'ward.D2'` explicitly |
| Same code, different dendrogram across R versions | R 3.1 renamed `'ward'` to `'ward.D'` and added `'ward.D2'` | Always specify `ward.D2` explicitly; never `'ward'` |
| Z-scored heatmap with extreme colors only in a few cells | Sparse matrix with zero-inflation | Filter low-expression rows; OR shift to robust scaling (`(x - median) / mad`) |
| Modules cluster by batch | Correlation distance picked up batch effect | Batch-correct upstream; verify via PCA |
| seaborn clustermap produces different clusters than R | seaborn `method='ward'` calls scipy.cluster.hierarchy.linkage which IS ward.D2-equivalent; difference is usually `metric` default | Set `metric='euclidean'` explicitly in both |
| OncoPrint mutual-exclusivity panel appears empty | `column_order` is being computed by clustering instead of preserved | Pass `column_order = ...` explicitly to oncoPrint |

**Operational rule:** a clustered heatmap is reproducible only when scaling, distance, linkage, color bounds, and (for OLO) the seriation method are all explicitly stated. Defaults differ across packages and across versions of the same package.

## Quantitative Thresholds

| Threshold | Value | Source |
|-----------|-------|--------|
| Robust color bound | 1st-99th percentile of |matrix| | Standard publication practice; suppresses single-outlier dominance |
| Raster trigger | >2000 rows or >2000 columns | ComplexHeatmap default `use_raster = TRUE` above 2000 |
| OLO practical limit | ~5000 rows | O(n^4) worst case; modern Bar-Joseph implementations faster |
| z-score symmetry | bounds around 0 | Z-scores are symmetric by construction |
| Minimum non-zero count to z-score | >=3 non-zero values per row | Below this sd is unreliable |
| Single-cell pseudobulk threshold | Downsample to <500 cells/group | Otherwise PDF rendering hangs |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| No PDF output from a script | Bare `Heatmap()` without `draw()` | Wrap in `draw()`; use `pdf()` / `dev.off()` explicitly |
| Color scale washed out | One outlier dominates | Clip to 1st-99th percentile; symmetric bounds for diverging |
| Time-course columns scrambled | `cluster_columns = TRUE` on ordered data | `cluster_columns = FALSE`; use `column_split` |
| pheatmap `gaps_col` ignored | Conflict with `cluster_cols = TRUE` | Disable clustering OR switch to ComplexHeatmap split |
| Dendrogram differs from a paper | Default `clustering_method` mismatch | Always specify `ward.D2`; never `ward` |
| seaborn standard_scale interpreted as z-score | Different rescaling functions | Use `z_score=0` for row z-scoring; `standard_scale` is min-max not z |
| Heatmap renders pixelated | Default raster_quality = 1 | Set `raster_quality = 5` for publication |
| Z-score blows up for some rows | Sparse rows; near-zero sd | Filter low-expression rows OR robust scale |
| Modules cluster by batch not biology | Batch effect not removed | limma::removeBatchEffect or ComBat upstream |

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Which clustering method?" | Explicit ward.D2 (Murtagh-Legendre 2014) on row-z-scored Euclidean distance. Alternatives evaluated in supplementary |
| "How were leaves ordered?" | Optimal Leaf Ordering via seriation::seriate (Bar-Joseph 2001); reduces visually-adjacent dissimilarity |
| "Why this color scale?" | Diverging palette symmetric around zero, bounds = 1st-99th percentile of |z|. Robust to outliers per standard practice |
| "Why z-score?" | Removes absolute level so co-regulated genes cluster regardless of magnitude. Raw values clustered separately (supplementary) |
| "Why row split by pathway?" | Pre-specified gene-set annotation (KEGG/Reactome) to verify clustering recovers known biology, not to bias it |
| "Reproducibility across R versions?" | `ward.D2` is stable across R 3.1+; `clustering_method = 'ward'` is not — never used |

## References

- Bar-Joseph Z, Gifford DK, Jaakkola TS. 2001. Fast optimal leaf ordering for hierarchical clustering. *Bioinformatics* 17(suppl 1):S22-S29. doi:10.1093/bioinformatics/17.suppl_1.S22
- Cerami E, Gao J, Dogrusoz U, et al. 2012. The cBio cancer genomics portal: an open platform for exploring multidimensional cancer genomics data. *Cancer Discov* 2(5):401-404.
- Crameri F, Shephard GE, Heron PJ. 2020. The misuse of colour in science communication. *Nat Commun* 11:5444.
- Gehlenborg N, Wong B. 2012. Points of view: Heat maps. *Nat Methods* 9(3):213.
- Gehlenborg N, Wong B. 2012. Points of view: Mapping quantitative data to color. *Nat Methods* 9(8):769.
- Gu Z, Eils R, Schlesner M. 2016. Complex heatmaps reveal patterns and correlations in multidimensional genomic data. *Bioinformatics* 32(18):2847-2849.
- Murtagh F, Legendre P. 2014. Ward's hierarchical agglomerative clustering method: which algorithms implement Ward's criterion? *J Classif* 31(3):274-295.
- Nuñez JR, Anderton CR, Renslow RS. 2018. Optimizing colormaps with consideration for color vision deficiency to enable accurate interpretation of scientific data. *PLOS ONE* 13(7):e0199239.

## Related Skills

- data-visualization/color-palettes - Sequential and diverging colormap selection
- data-visualization/oncoprint-mutation-matrices - Mutation-matrix heatmap (ComplexHeatmap oncoPrint)
- data-visualization/multipanel-figures - Combine heatmaps into journal layouts
- data-visualization/dimensionality-reduction-plots - PCA / UMAP as alternative views of the same matrix
- differential-expression/de-visualization - Heatmap of top DE genes
- single-cell/markers-annotation - Single-cell dotplot / matrixplot as scRNA alternatives
