# Heatmaps and Clustering Usage Guide

## Overview

Clustered heatmaps visualize matrix data (expression, methylation, etc.) with hierarchical clustering to reveal patterns across samples and features.

## Tool Selection

| Tool | Language | Best For |
|------|----------|----------|
| pheatmap | R | Quick, publication-ready heatmaps |
| ComplexHeatmap | R | Complex annotations, multiple heatmaps |
| seaborn clustermap | Python | Python workflows, simple annotations |

## Quick Start Prompts

- "Create a heatmap of my top DE genes with sample annotations"
- "Cluster my expression matrix and extract cluster assignments"
- "Make a heatmap with row and column annotations"
- "Split my heatmap by gene pathway"

## Requirements

```r
# R
install.packages('pheatmap')
BiocManager::install('ComplexHeatmap')
```

```bash
# Python
pip install seaborn matplotlib scipy
```

## Key Considerations

- **Scale data** - Usually z-score rows for expression
- **Distance metric** - Correlation for expression, euclidean for abundance
- **Linkage method** - ward.D2 often works well
- **Color palette** - Diverging for centered data (RdBu), sequential otherwise

## Related Skills

- **data-visualization/color-palettes** - Color selection
- **differential-expression/de-visualization** - DE-specific plots
