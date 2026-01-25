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

- Always use variance-stabilized counts (vst) for PCA and heatmaps
- Scale heatmap rows (z-score) for comparable gene patterns
- Check p-value histogram for analysis quality
- Use colorblind-friendly palettes for publications
- Save vector formats (PDF) for publications, raster (PNG) for presentations
