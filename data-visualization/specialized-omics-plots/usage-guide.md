# Specialized Omics Plots - Usage Guide

## Overview
Standard omics visualizations for differential expression, dimensionality reduction, and enrichment results.

## Prerequisites
```r
# R
install.packages(c('ggplot2', 'ggrepel', 'ggpubr', 'corrplot'))
BiocManager::install('survminer')
```

```bash
# Python
pip install matplotlib seaborn scikit-learn scanpy
```

## Quick Start
Tell your AI agent what you want to do:
- "Create a volcano plot with labeled significant genes"
- "Make a PCA colored by condition with ellipses"
- "Plot enrichment results as a dotplot"

## Example Prompts
### Differential Expression
> "Create a volcano plot highlighting genes with FDR < 0.05 and |log2FC| > 1"

> "Make an MA plot with significant genes colored"

### Dimensionality Reduction
> "Create a PCA plot colored by treatment with 95% confidence ellipses"

> "Make a UMAP of my single-cell data colored by cluster"

### Enrichment Visualization
> "Create a dotplot of my GO enrichment results"

> "Make an enrichment barplot sorted by significance"

### Statistical Plots
> "Compare expression across groups with boxplots and significance stars"

> "Create a correlation heatmap of my samples"

## What the Agent Will Do
1. Prepare data in appropriate format
2. Create the visualization with proper aesthetics
3. Add labels, annotations, and statistical marks
4. Apply publication-ready theme
5. Export at high resolution

## Common Plot Types

| Plot | Use Case |
|------|----------|
| Volcano | DE significance vs fold change |
| MA | Expression level vs fold change |
| PCA | Sample relationships, batch effects |
| Dotplot | Enrichment results |
| UMAP/tSNE | Single-cell clustering |
| Heatmap | Gene expression patterns |

## Tips
- Volcano: Use diverging colors (up=red, down=blue)
- PCA: Include variance explained in axis labels
- Dotplot: Sort by significance, size by count
- Always add appropriate statistical annotations

## Related Skills
- **data-visualization/ggplot2-fundamentals** - Base syntax
- **data-visualization/color-palettes** - Colors
