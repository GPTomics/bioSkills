# KEGG Pathway Enrichment - Usage Guide

## Overview
KEGG (Kyoto Encyclopedia of Genes and Genomes) pathway enrichment tests whether your genes are enriched in specific biological pathways including signaling, metabolism, and disease pathways.

## Prerequisites
```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install(c('clusterProfiler', 'org.Hs.eg.db'))
```

## Quick Start
Tell your AI agent what you want to do:
- "Run KEGG pathway enrichment on my significant genes"
- "Find which metabolic pathways are affected in my experiment"
- "Identify signaling pathways enriched in my gene list"

## Example Prompts
### Basic Enrichment
> "Run KEGG pathway enrichment on my differentially expressed genes from deseq2_results.csv"

> "Find enriched KEGG pathways in my significant genes using human annotation"

### Organism-Specific
> "Run KEGG enrichment on my mouse gene list using mmu organism code"

> "What organism code should I use for zebrafish KEGG analysis?"

### Visualization
> "Run KEGG enrichment and show a dotplot of the top 20 pathways"

> "Open the most enriched KEGG pathway in the browser"

### KEGG Modules
> "Run KEGG module enrichment instead of full pathway enrichment"

## What the Agent Will Do
1. Load DE results and extract significant genes
2. Convert gene IDs to Entrez format (required for KEGG)
3. Run enrichKEGG() with appropriate organism code
4. Convert results to readable gene symbols with setReadable()
5. Generate visualizations and export results

## Common Organism Codes

| Code | Organism |
|------|----------|
| hsa | Human |
| mmu | Mouse |
| rno | Rat |
| dre | Zebrafish |
| dme | Drosophila |
| cel | C. elegans |
| sce | S. cerevisiae |
| ath | Arabidopsis |

## Understanding Results

| Column | Description |
|--------|-------------|
| ID | KEGG pathway ID (hsa04110) |
| Description | Pathway name |
| GeneRatio | Genes in pathway / Total query |
| BgRatio | Background genes in pathway / Total |
| pvalue | Raw p-value |
| p.adjust | FDR-adjusted p-value |
| geneID | Genes in the pathway |
| Count | Number of genes |

## Tips
- KEGG requires Entrez gene IDs - always convert from symbols first
- Use search_kegg_organism() to find the correct organism code
- Use setReadable() to convert Entrez IDs back to symbols in results
- KEGG queries require internet connection; use use_internal_data = TRUE for cached data
- See enrichment-visualization skill for plotting (dotplot, barplot, browseKEGG)
- Try enrichMKEGG() for KEGG module analysis (smaller functional units)
