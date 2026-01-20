# Enrichment Visualization - Usage Guide

## Overview
The enrichplot package provides visualization functions for clusterProfiler results, including dot plots, bar plots, networks, and GSEA-specific plots.

## Prerequisites
```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install(c('clusterProfiler', 'enrichplot'))
```

## Quick Start
Tell your AI agent what you want to do:
- "Create a dotplot of my GO enrichment results"
- "Make a gene-concept network from my KEGG enrichment"
- "Show a GSEA running score plot for my top pathway"

## Example Prompts
### Dot Plots
> "Create a dotplot showing the top 20 enriched GO terms"

> "Make a dotplot of my KEGG results with adjusted font size for long pathway names"

### Network Plots
> "Create a gene-concept network colored by fold change"

> "Generate an enrichment map clustering similar GO terms together"

### GSEA Plots
> "Show a GSEA running score plot for the top 3 enriched pathways"

> "Create a ridge plot showing fold change distributions for each gene set"

### Customization
> "Save the enrichment dotplot as a PDF at publication quality"

> "Change the color scale to viridis on my enrichment plot"

### Comparison Plots
> "Create a dotplot comparing enrichment results between up and down regulated genes"

## What the Agent Will Do
1. Take enrichment results from GO, KEGG, or GSEA analysis
2. Select appropriate plot type based on request
3. For emapplot, compute pairwise term similarity first
4. Generate publication-quality figure with appropriate sizing
5. Save to PDF or PNG with specified dimensions

## Plot Types Quick Reference

| Plot Type | Function | Best For |
|-----------|----------|----------|
| Dot plot | dotplot() | Overview of top terms |
| Bar plot | barplot() | Simple count/ratio display |
| Network | cnetplot() | Gene-concept relationships |
| Map | emapplot() | Term similarity clusters |
| Tree | treeplot() | Hierarchical term grouping |
| Upset | upsetplot() | Overlapping genes |
| GSEA | gseaplot2() | Running enrichment score |
| Ridge | ridgeplot() | Fold change distribution |
| Heatmap | heatplot() | Gene-concept matrix |

## Choosing a Visualization

### For Over-Representation Results
- Starting point: dotplot() - best overview
- Show relationships: cnetplot() - genes and terms
- Cluster terms: emapplot() - similar terms grouped
- Compare groups: dotplot() with compareCluster result

### For GSEA Results
- Starting point: ridgeplot() - all gene sets
- Detailed view: gseaplot2() - running score for specific sets
- Overview: dotplot() - works for GSEA too

## Tips
- dotplot() is the most common choice for publications - start there
- For emapplot(), always run pairwise_termsim() first on your enrichment result
- Use showCategory parameter to control number of terms displayed
- For long term names, increase plot width and reduce font size
- cnetplot() with circular = TRUE helps with crowded networks
- gseaplot2() accepts multiple geneSetID values to compare pathways
- Use ggsave() for better control over output dimensions and resolution
