# GO Enrichment Usage Guide

Gene Ontology over-representation analysis tests whether specific GO terms appear more frequently in your gene list than expected by chance.

## When to Use

- You have a list of significant genes from differential expression
- You want to understand biological processes, molecular functions, or cellular components
- You need statistical evidence that specific functions are enriched

## Prerequisites

```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install(c('clusterProfiler', 'org.Hs.eg.db'))
```

## Basic Workflow

### 1. Prepare Your Gene List

Gene list should be a character vector of gene IDs:

```r
# From DE results
de_results <- read.csv('deseq2_results.csv')
sig_genes <- de_results$gene_id[de_results$padj < 0.05 & abs(de_results$log2FoldChange) > 1]
```

### 2. Convert Gene IDs (if needed)

```r
library(clusterProfiler)
library(org.Hs.eg.db)

# Convert symbols to Entrez IDs
gene_ids <- bitr(sig_genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
gene_list <- gene_ids$ENTREZID
```

### 3. Run Enrichment

```r
ego <- enrichGO(
    gene = gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = 'ENTREZID',
    ont = 'BP',  # Biological Process
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
)
```

### 4. View Results

```r
# Summary
head(ego)

# Full data frame
results <- as.data.frame(ego)
```

## Understanding Results

| Column | Description |
|--------|-------------|
| ID | GO term ID (GO:XXXXXXX) |
| Description | GO term name |
| GeneRatio | Genes in term / Total query genes |
| BgRatio | Background genes in term / Total background |
| pvalue | Raw p-value (hypergeometric test) |
| p.adjust | Adjusted p-value (FDR) |
| qvalue | Q-value |
| geneID | Genes in the term |
| Count | Number of genes |

## Three Ontologies

| Ontology | Code | Description |
|----------|------|-------------|
| Biological Process | BP | What the genes do |
| Molecular Function | MF | Biochemical activity |
| Cellular Component | CC | Where in cell |

## Tips

### Use a Background Universe

Without specifying a background, clusterProfiler uses all genes in the OrgDb. For better results, use your expressed genes:

```r
all_expressed <- de_results$gene_id  # All genes tested
universe <- bitr(all_expressed, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)

ego <- enrichGO(gene = gene_list, universe = universe$ENTREZID, ...)
```

### Simplify Redundant Terms

GO terms are hierarchical, causing redundant enrichments:

```r
ego_simple <- simplify(ego, cutoff = 0.7)
```

### Separate Up/Down Genes

Run enrichment separately on up- and down-regulated genes:

```r
up_genes <- de_results$gene_id[de_results$padj < 0.05 & de_results$log2FoldChange > 1]
down_genes <- de_results$gene_id[de_results$padj < 0.05 & de_results$log2FoldChange < -1]

ego_up <- enrichGO(gene = up_genes_entrez, ...)
ego_down <- enrichGO(gene = down_genes_entrez, ...)
```

## Common Issues

**No enriched terms found:**
- Check gene ID conversion (many IDs may fail to convert)
- Loosen thresholds (increase pvalueCutoff)
- Check organism database matches your data

**Too many terms:**
- Use `simplify()` to reduce redundancy
- Lower pvalueCutoff or qvalueCutoff
- Increase minGSSize

## Visualization

See the enrichment-visualization skill for plotting options:
- `dotplot(ego)` - Dot plot of top terms
- `barplot(ego)` - Bar plot of enrichment
- `cnetplot(ego)` - Gene-concept network
- `emapplot(ego)` - Enrichment map
