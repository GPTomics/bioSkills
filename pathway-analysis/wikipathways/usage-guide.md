# WikiPathways Enrichment - Usage Guide

## Overview

WikiPathways is an open, collaborative platform for biological pathways. Unlike KEGG (commercial license) or Reactome (limited species), WikiPathways is community-curated with CC0 license and supports 30+ species. The clusterProfiler package provides enrichWP() and gseWP() functions for direct analysis.

## When to Use This Skill

- You want open-source pathway annotations (CC0 license)
- You're working with species not in Reactome
- You want disease-specific or drug-related pathways
- You need community-contributed specialized pathways
- You want to complement KEGG/Reactome results

## Installation

```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install(c('clusterProfiler', 'rWikiPathways', 'enrichplot'))
BiocManager::install('org.Hs.eg.db')  # For human
```

## Basic Workflow

### 1. Prepare Gene List

```r
library(clusterProfiler)
library(org.Hs.eg.db)

de_results <- read.csv('deseq2_results.csv')
sig_genes <- de_results[de_results$padj < 0.05 & abs(de_results$log2FoldChange) > 1, ]

gene_ids <- bitr(sig_genes$gene_symbol, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
entrez_ids <- gene_ids$ENTREZID
```

### 2. Run Enrichment

```r
wp_result <- enrichWP(
    gene = entrez_ids,
    organism = 'Homo sapiens',
    pvalueCutoff = 0.05
)

# Make readable
wp_readable <- setReadable(wp_result, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
head(as.data.frame(wp_readable))
```

### 3. Visualize

```r
library(enrichplot)

dotplot(wp_readable, showCategory = 15)
```

## Exploring WikiPathways

```r
library(rWikiPathways)

# Available organisms
listOrganisms()

# Pathways for organism
pathways <- listPathways('Homo sapiens')
head(pathways)

# Search pathways
searchPathways('cancer', 'Homo sapiens')
searchPathways('insulin signaling', 'Homo sapiens')

# Pathway details
getPathwayInfo('WP554')
```

## Common Issues

### Organism Name

Use the full scientific name exactly as listed:
```r
library(rWikiPathways)
listOrganisms()  # Check exact spelling
```

### Internet Connection

WikiPathways requires internet access to download current pathway data.

### Few Results

WikiPathways has fewer pathways than KEGG. Consider:
- Using multiple databases (WikiPathways + KEGG + Reactome)
- Relaxing p-value cutoff
- Checking if pathways exist for your research area

## Output Format

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

## Resources

- [WikiPathways Website](https://www.wikipathways.org/)
- [rWikiPathways Bioconductor](https://bioconductor.org/packages/rWikiPathways/)
- [WikiPathways API](https://webservice.wikipathways.org/)
