# Heatmaps and Clustering - Usage Guide

## Overview
Clustered heatmaps visualize matrix data (expression, methylation, etc.) with hierarchical clustering to reveal patterns across samples and features.

## Prerequisites
```r
# R
install.packages('pheatmap')
BiocManager::install('ComplexHeatmap')
```

```bash
# Python
pip install seaborn matplotlib scipy
```

## Quick Start
Tell your AI agent what you want to do:
- "Create a heatmap of my top DE genes with sample annotations"
- "Cluster my expression matrix and extract cluster assignments"
- "Make a heatmap with row and column annotations"

## Example Prompts
### Basic Heatmaps
> "Create a heatmap of my expression matrix with clustering"

> "Make a quick heatmap of my top 50 variable genes"

### Annotations
> "Add sample annotations showing treatment and batch"

> "Color row annotations by gene pathway"

### Complex Layouts
> "Split my heatmap by gene pathway"

> "Create side-by-side heatmaps for two conditions"

## What the Agent Will Do
1. Prepare and scale the data matrix (typically z-score rows)
2. Set up annotation data frames for rows/columns
3. Choose appropriate distance metric and clustering method
4. Configure color palette and legend
5. Generate and export the heatmap

## Tool Selection

| Tool | Language | Best For |
|------|----------|----------|
| pheatmap | R | Quick, publication-ready heatmaps |
| ComplexHeatmap | R | Complex annotations, multiple heatmaps |
| seaborn clustermap | Python | Python workflows, simple annotations |

## Tips
- Scale data (z-score rows) for expression data
- Use correlation distance for expression, euclidean for abundance
- ward.D2 linkage often works well for biological data
- Use diverging palettes (RdBu) for centered data, sequential otherwise
- Limit to 50-100 genes for readability

## Related Skills
- **data-visualization/color-palettes** - Color selection
- **differential-expression/de-visualization** - DE-specific plots
