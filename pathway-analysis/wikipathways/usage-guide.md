# WikiPathways Enrichment - Usage Guide

## Overview
WikiPathways is an open, collaborative platform for biological pathways with CC0 license, community curation, and support for 30+ species including many not covered by KEGG or Reactome.

## Prerequisites
```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install(c('clusterProfiler', 'rWikiPathways', 'enrichplot', 'org.Hs.eg.db'))
```

## Quick Start
Tell your AI agent what you want to do:
- "Run WikiPathways enrichment on my significant genes"
- "Find disease-specific pathways in my gene list"
- "Use open-source pathway database for enrichment"

## Example Prompts
### Basic Enrichment
> "Run WikiPathways enrichment on my differentially expressed genes from deseq2_results.csv"

> "Perform WikiPathways enrichment on my significant genes for human"

### Organism-Specific
> "Run WikiPathways enrichment for mouse genes"

> "What organisms are available in WikiPathways?"

### Exploring Pathways
> "Search WikiPathways for cancer-related pathways"

> "Find drug metabolism pathways in WikiPathways"

### Combining Databases
> "Run enrichment using WikiPathways, KEGG, and Reactome and compare results"

### GSEA
> "Run WikiPathways GSEA on my ranked gene list"

## What the Agent Will Do
1. Load DE results and extract significant genes
2. Convert gene IDs to Entrez format
3. Run enrichWP() with appropriate organism name
4. Convert results to readable gene symbols with setReadable()
5. Generate visualizations and export results

## Exploring Available Pathways
```r
library(rWikiPathways)

# List available organisms
listOrganisms()

# Search for pathways
searchPathways('cancer', 'Homo sapiens')
searchPathways('insulin signaling', 'Homo sapiens')

# Get pathway details
getPathwayInfo('WP554')
```

## WikiPathways vs Other Databases

| Feature | WikiPathways | KEGG | Reactome |
|---------|--------------|------|----------|
| License | CC0 (open) | Commercial | Open |
| Curation | Community | Expert | Peer-reviewed |
| Species | 30+ | 4000+ | 7 |
| Focus | Disease/drug pathways | Metabolic/signaling | Reactions |

## Understanding Results

| Column | Description |
|--------|-------------|
| ID | WikiPathways ID (WP####) |
| Description | Pathway name |
| GeneRatio | Genes in list / genes in pathway |
| BgRatio | Pathway genes / total genes |
| pvalue | Raw p-value |
| p.adjust | Adjusted p-value |
| geneID | Genes in pathway |
| Count | Number of genes |

## Tips
- Use the exact scientific name from listOrganisms() for the organism parameter
- WikiPathways requires internet access to download current pathway data
- WikiPathways has fewer total pathways than KEGG - combine with other databases
- Community pathways include specialized disease and drug-related pathways
- Use setReadable() to convert Entrez IDs to gene symbols in results
- See enrichment-visualization skill for dotplot() and other plots
