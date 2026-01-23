---
name: bio-pathway-go-enrichment
description: Gene Ontology over-representation analysis using clusterProfiler enrichGO. Use when identifying biological functions enriched in a gene list from differential expression or other analyses. Supports all three ontologies (BP, MF, CC), multiple ID types, and customizable statistical thresholds.
tool_type: r
primary_tool: clusterProfiler
---

# GO Over-Representation Analysis

## Core Pattern

```r
library(clusterProfiler)
library(org.Hs.eg.db)  # Human - change for other organisms

ego <- enrichGO(
    gene = gene_list,           # Character vector of gene IDs
    OrgDb = org.Hs.eg.db,       # Organism annotation database
    keyType = 'ENTREZID',       # ID type: ENSEMBL, SYMBOL, ENTREZID, etc.
    ont = 'BP',                 # BP, MF, CC, or ALL
    pAdjustMethod = 'BH',       # p-value adjustment method
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
)
```

## Prepare Gene List from DE Results

```r
library(dplyr)

de_results <- read.csv('de_results.csv')

sig_genes <- de_results %>%
    filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
    pull(gene_id)

# If using gene symbols, convert to Entrez IDs
gene_ids <- bitr(sig_genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
gene_list <- gene_ids$ENTREZID
```

## ID Conversion with bitr

```r
# Check available key types
keytypes(org.Hs.eg.db)

# Convert between ID types
converted <- bitr(genes, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)

# Multiple output types
converted <- bitr(genes, fromType = 'SYMBOL', toType = c('ENTREZID', 'ENSEMBL'), OrgDb = org.Hs.eg.db)
```

## With Background Universe

```r
# Use all expressed genes as background (recommended)
all_genes <- de_results$gene_id
universe_ids <- bitr(all_genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)

ego <- enrichGO(
    gene = gene_list,
    universe = universe_ids$ENTREZID,  # Background gene set
    OrgDb = org.Hs.eg.db,
    keyType = 'ENTREZID',
    ont = 'BP',
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05
)
```

## All Three Ontologies

```r
# Run all ontologies at once
ego_all <- enrichGO(
    gene = gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = 'ENTREZID',
    ont = 'ALL',  # BP, MF, and CC combined
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05
)

# Results include ONTOLOGY column
head(as.data.frame(ego_all))
```

## Make Results Readable

```r
# Convert Entrez IDs to gene symbols in results
ego_readable <- setReadable(ego, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')

# Or use readable = TRUE directly (only works with ENTREZID input)
ego <- enrichGO(
    gene = gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = 'ENTREZID',
    ont = 'BP',
    readable = TRUE  # Converts to symbols
)
```

## Extract and Export Results

```r
# View top results
head(ego)

# Convert to data frame
results_df <- as.data.frame(ego)

# Key columns: ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count

# Export to CSV
write.csv(results_df, 'go_enrichment_results.csv', row.names = FALSE)

# Filter for specific criteria
sig_terms <- results_df[results_df$p.adjust < 0.01 & results_df$Count >= 5, ]
```

## Simplify Redundant Terms

```r
# Remove redundant GO terms (keeps representative terms)
ego_simplified <- simplify(ego, cutoff = 0.7, by = 'p.adjust', select_fun = min)
```

## Different Organisms

```r
# Mouse
library(org.Mm.eg.db)
ego_mouse <- enrichGO(gene = genes, OrgDb = org.Mm.eg.db, ont = 'BP')

# Zebrafish
library(org.Dr.eg.db)
ego_zfish <- enrichGO(gene = genes, OrgDb = org.Dr.eg.db, ont = 'BP')

# Yeast
library(org.Sc.sgd.db)
ego_yeast <- enrichGO(gene = genes, OrgDb = org.Sc.sgd.db, ont = 'BP', keyType = 'ORF')
```

## Group GO Terms by Ancestor

```r
# Classify genes by GO slim categories
ggo <- groupGO(
    gene = gene_list,
    OrgDb = org.Hs.eg.db,
    ont = 'BP',
    level = 3,  # GO hierarchy level
    readable = TRUE
)
```

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| gene | required | Vector of gene IDs |
| OrgDb | required | Organism database |
| keyType | ENTREZID | Input ID type |
| ont | BP | BP, MF, CC, or ALL |
| pvalueCutoff | 0.05 | P-value threshold |
| qvalueCutoff | 0.2 | Q-value (FDR) threshold |
| pAdjustMethod | BH | BH, bonferroni, etc. |
| universe | NULL | Background genes |
| minGSSize | 10 | Min genes per term |
| maxGSSize | 500 | Max genes per term |
| readable | FALSE | Convert to symbols |

## Related Skills

- kegg-pathways - KEGG pathway enrichment
- gsea - Gene Set Enrichment Analysis for GO
- enrichment-visualization - Visualize enrichment results
- differential-expression - Generate input gene lists
