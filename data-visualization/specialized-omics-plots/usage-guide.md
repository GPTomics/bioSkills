# Specialized Omics Plots Usage Guide

## Overview

Standard omics visualizations for differential expression, dimensionality reduction, and enrichment results.

## Common Plot Types

| Plot | Use Case |
|------|----------|
| Volcano | DE significance vs fold change |
| MA | Expression level vs fold change |
| PCA | Sample relationships, batch effects |
| Dotplot | Enrichment results |
| UMAP/tSNE | Single-cell clustering |
| Heatmap | Gene expression patterns |

## Quick Start Prompts

- "Create a volcano plot with labeled significant genes"
- "Make a PCA colored by condition with ellipses"
- "Plot enrichment results as a dotplot"
- "Compare expression across groups with boxplots"

## Requirements

```r
# R
install.packages(c('ggplot2', 'ggrepel', 'ggpubr', 'corrplot'))
BiocManager::install('survminer')
```

```bash
# Python
pip install matplotlib seaborn scikit-learn scanpy
```

## Key Principles

- **Volcano**: Use diverging colors (up=red, down=blue)
- **PCA**: Include variance explained in axis labels
- **Dotplot**: Sort by significance, size by count
- **All**: Add appropriate statistical annotations

## Related Skills

- **data-visualization/ggplot2-fundamentals** - Base syntax
- **data-visualization/color-palettes** - Colors
