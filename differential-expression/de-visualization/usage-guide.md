# DE Visualization - Usage Guide

## Overview

This skill covers creating publication-quality visualizations for differential expression results, including MA plots, volcano plots, PCA plots, and heatmaps. Works with both DESeq2 and edgeR output.

## Prerequisites

```r
install.packages(c('ggplot2', 'pheatmap', 'RColorBrewer', 'ggrepel'))
BiocManager::install('EnhancedVolcano')  # Optional
```

## Quick Start

Tell your AI agent what you want to do:
- "Create a volcano plot from my DESeq2 results"
- "Make a heatmap of the top 50 differentially expressed genes"
- "Generate a PCA plot colored by treatment group"

## Example Prompts

### MA and Volcano Plots
> "Create an MA plot highlighting significant genes"

> "Make a volcano plot with gene labels for top hits"

> "Generate an EnhancedVolcano plot for my results"

### PCA and Sample Clustering
> "Show PCA of my samples colored by condition and shaped by batch"

> "Create a sample distance heatmap"

> "Plot MDS for my edgeR data"

### Heatmaps
> "Make a heatmap of significant genes clustered by expression"

> "Create a heatmap for my genes of interest"

> "Show expression patterns across samples"

### Individual Genes
> "Plot counts for gene X across conditions"

> "Show expression of my candidate genes"

## What the Agent Will Do

1. Extract and format results from DESeq2/edgeR
2. Apply appropriate transformations (vst, log)
3. Create publication-quality figures
4. Add annotations and labels
5. Save in requested format (PDF, PNG)

## Common Plot Types

| Plot | Shows | Use For |
|------|-------|---------|
| MA plot | LFC vs expression | QC, global view |
| Volcano | LFC vs significance | Identifying top genes |
| PCA | Sample relationships | Batch effects, outliers |
| Heatmap | Expression patterns | Gene clusters, validation |

## Tips

- Always use variance-stabilized counts (vst or rlog) for PCA and heatmaps; raw counts make PC1 = library size
- `vst()` function default is `blind=TRUE`; the current DESeq2 vignette recommends `blind=FALSE` for downstream visualization after the model is fit -- pass `blind=` explicitly
- Scale heatmap rows (z-score) for comparable gene patterns; for QC heatmaps showing batch shifts, use `scale='none'` instead
- Check p-value histogram for analysis quality (uniform + spike near 0 is correct; U-shape means batch effects)
- Use shrunken LFCs for volcano plot x-axis and un-shrunken p-values for y-axis; shrinkage compresses the cloud but does NOT recompute p-values, so the same genes remain significant
- Always set `max.overlaps = Inf` (or `options(ggrepel.max.overlaps = Inf)`) when labeling more than ~10 genes -- ggrepel's default silently drops labels
- For outlier-robust top-variable-gene selection, use `matrixStats::rowMads(assay(vsd))` instead of `rowVars` -- one outlier sample can inflate rowVars massively
- Use colorblind-friendly palettes for publications
- Save vector formats (PDF) for publications, raster (PNG) for presentations
- MA plot cloud should be symmetric around LFC=0; asymmetry suggests normalization failure (TMM/RLE assumption violated)

## Related Skills

- deseq2-basics - Generate DESeq2 results for visualization
- edger-basics - Generate edgeR results for visualization
- de-results - Filter genes before visualization
- data-visualization/volcano-and-ma-plots - Custom ggplot2/matplotlib volcano + MA with LFC shrinkage
- data-visualization/dimensionality-reduction-plots - PCA / UMAP / t-SNE / PHATE
- data-visualization/heatmaps-clustering - Advanced heatmap customization
