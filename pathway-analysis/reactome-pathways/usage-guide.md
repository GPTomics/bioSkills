# Reactome Pathway Enrichment - Usage Guide

## Overview
ReactomePA provides pathway enrichment analysis using the Reactome database, a curated peer-reviewed knowledgebase of biological pathways with deep hierarchy and reaction-level detail.

## Prerequisites
```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install(c('ReactomePA', 'org.Hs.eg.db', 'enrichplot'))
```

## Quick Start
Tell your AI agent what you want to do:
- "Run Reactome pathway enrichment on my significant genes"
- "Find peer-reviewed curated pathways enriched in my gene list"
- "Use Reactome instead of KEGG for pathway analysis"

## Example Prompts
### Basic Enrichment
> "Run Reactome pathway enrichment on my differentially expressed genes from deseq2_results.csv"

> "Perform Reactome enrichment on my significant genes and show the top 20 pathways"

### Visualization
> "Create a dotplot of my Reactome enrichment results"

> "Open the most enriched Reactome pathway in the browser"

### Comparison with KEGG
> "Run both Reactome and KEGG enrichment on my genes and compare results"

### GSEA with Reactome
> "Run Reactome GSEA on my ranked gene list"

## What the Agent Will Do
1. Load DE results and extract significant genes
2. Convert gene IDs to Entrez format (required for ReactomePA)
3. Run enrichPathway() with readable = TRUE for gene symbols
4. Generate visualizations and export results
5. Optionally view pathways in Reactome browser

## Reactome vs KEGG

| Feature | Reactome | KEGG |
|---------|----------|------|
| Curation | Peer-reviewed | Expert curated |
| Access | Open source | Requires license for commercial |
| Hierarchy | Deep pathway hierarchy | Flat pathway list |
| Detail | Reaction-level detail | Pathway-level |
| Updates | Quarterly | Ongoing |
| Organisms | 7 species | 4000+ species |

## Understanding Results

| Column | Description |
|--------|-------------|
| ID | Reactome pathway ID (R-HSA-XXXXX) |
| Description | Pathway name |
| GeneRatio | Genes in list / genes in pathway |
| BgRatio | Pathway genes / total genes |
| pvalue | Raw p-value |
| p.adjust | Adjusted p-value (BH) |
| qvalue | Q-value |
| geneID | Genes in pathway |
| Count | Number of genes |

## Tips
- Reactome requires Entrez gene IDs - convert from symbols first with bitr()
- Use readable = TRUE in enrichPathway() to get gene symbols in output
- viewPathway() opens interactive pathway visualization in browser
- Reactome has fewer species than KEGG but more detailed pathway annotation
- Combine Reactome with KEGG/WikiPathways for comprehensive pathway coverage
- See enrichment-visualization skill for dotplot() and other plots
