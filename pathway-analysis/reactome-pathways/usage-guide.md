# Reactome Pathway Enrichment - Usage Guide

## Overview

ReactomePA provides pathway enrichment analysis using the Reactome database, a curated peer-reviewed knowledgebase of biological pathways. It offers both over-representation analysis (ORA) and Gene Set Enrichment Analysis (GSEA) for Reactome pathways.

## When to Use This Skill

- You have a list of differentially expressed genes and want pathway enrichment
- You prefer curated, peer-reviewed pathway annotations
- You want detailed pathway hierarchy information
- You need to visualize pathways in the Reactome browser
- You want an alternative to KEGG pathways

## Reactome vs KEGG

| Feature | Reactome | KEGG |
|---------|----------|------|
| Curation | Peer-reviewed | Expert curated |
| Access | Open source | Requires license for commercial |
| Hierarchy | Deep pathway hierarchy | Flat pathway list |
| Detail | Reaction-level detail | Pathway-level |
| Updates | Quarterly | Ongoing |
| Organisms | 7 species | 4000+ species |

## Installation

```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install('ReactomePA')

# Also need organism database
BiocManager::install('org.Hs.eg.db')  # Human
```

## Basic Workflow

### 1. Prepare Gene List

```r
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)

# From differential expression results
de_results <- read.csv('deseq2_results.csv')

# Get significant genes
sig_genes <- de_results[de_results$padj < 0.05 & abs(de_results$log2FoldChange) > 1, ]

# Convert to Entrez IDs (required for ReactomePA)
gene_ids <- bitr(sig_genes$gene_symbol, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
entrez_ids <- gene_ids$ENTREZID
```

### 2. Run Enrichment

```r
pathway_result <- enrichPathway(
    gene = entrez_ids,
    organism = 'human',
    pvalueCutoff = 0.05,
    readable = TRUE
)

# View results
head(as.data.frame(pathway_result))
```

### 3. Visualize

```r
library(enrichplot)

# Dot plot
dotplot(pathway_result, showCategory = 20)

# View specific pathway in browser
viewPathway(pathway_result@result$ID[1], organism = 'human')
```

## Common Issues

### No Results

- Check that genes are Entrez IDs (not symbols or Ensembl)
- Try relaxing p-value cutoff
- Ensure organism parameter matches your data
- Check that genes are in Reactome database

### ID Conversion Failures

```r
# Check conversion success rate
converted <- bitr(genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
message(sprintf('Converted %d of %d genes', nrow(converted), length(genes)))
```

## Output Format

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

## Resources

- [Reactome Website](https://reactome.org/)
- [ReactomePA Bioconductor](https://bioconductor.org/packages/ReactomePA/)
- [ReactomePA Documentation](https://yulab-smu.top/biomedical-knowledge-mining-book/reactomepa.html)
