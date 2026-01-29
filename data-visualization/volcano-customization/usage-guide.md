# Volcano Plot Customization - Usage Guide

## Overview

Volcano plots display statistical significance (-log10 p-value) vs effect size (log2 fold change) for differential expression or association results, with customizable thresholds and gene labels.

## Prerequisites

```r
# R
install.packages(c('ggplot2', 'ggrepel', 'dplyr'))
BiocManager::install('EnhancedVolcano')
```

```bash
# Python
pip install matplotlib seaborn adjustText pandas numpy statsmodels
```

## Quick Start

Tell your AI agent:
- "Create a volcano plot from my DE results"
- "Label the top 20 significant genes on my volcano plot"
- "Highlight specific genes of interest on the volcano"
- "Adjust thresholds to FC > 1.5 and padj < 0.01"

## Example Prompts

### Basic Volcano Plots

> "Create a volcano plot from my DESeq2 results"

> "Make a volcano showing upregulated genes in red, downregulated in blue"

### Gene Labeling

> "Label the top 10 most significant genes"

> "Highlight and label TP53, BRCA1, and MYC on the volcano"

> "Only label genes with FC > 2 and padj < 0.001"

### Customization

> "Change the fold change threshold to 1.5"

> "Add dashed lines at the significance thresholds"

> "Use a colorblind-friendly palette for the volcano"

## What the Agent Will Do

1. Load DE results and identify required columns (log2FC, pvalue/padj, gene)
2. Calculate significance categories based on thresholds
3. Create base volcano plot with appropriate colors
4. Add threshold indicator lines
5. Label requested genes using text repulsion
6. Export in requested format and resolution

## Tool Selection

| Tool | Language | Best For |
|------|----------|----------|
| ggplot2 + ggrepel | R | Full customization, publication figures |
| EnhancedVolcano | R | Quick, feature-rich volcanos |
| matplotlib + adjustText | Python | Python workflows |

## Tips

- **Use padj (adjusted p-value)** for coloring to account for multiple testing
- **ggrepel/adjustText** prevent label overlaps automatically
- **Limit labels** to top 15-25 genes for readability
- **FC threshold 1** = 2-fold change; 0.58 = 1.5-fold change
- **Consider asymmetric thresholds** if biology suggests one direction
- **Box labels** (EnhancedVolcano) improve readability with many labels
