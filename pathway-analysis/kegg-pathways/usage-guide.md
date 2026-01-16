# KEGG Pathway Enrichment Usage Guide

KEGG (Kyoto Encyclopedia of Genes and Genomes) pathway enrichment tests whether your genes are enriched in specific biological pathways.

## When to Use

- You have a list of significant genes
- You want to identify affected pathways (signaling, metabolism, disease)
- You need to interpret results in terms of pathway biology

## Prerequisites

```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install('clusterProfiler')

# For converting gene IDs:
BiocManager::install('org.Hs.eg.db')  # Human
```

## Basic Workflow

### 1. Prepare Gene List

KEGG requires Entrez gene IDs:

```r
library(clusterProfiler)
library(org.Hs.eg.db)

de_results <- read.csv('de_results.csv')
sig_genes <- de_results$gene_id[de_results$padj < 0.05 & abs(de_results$log2FoldChange) > 1]

# Convert to Entrez IDs
gene_ids <- bitr(sig_genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
gene_list <- gene_ids$ENTREZID
```

### 2. Run Enrichment

```r
kk <- enrichKEGG(
    gene = gene_list,
    organism = 'hsa',  # Human
    pvalueCutoff = 0.05
)
```

### 3. Make Readable

```r
# Convert Entrez IDs to symbols in results
kk_readable <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')

# View results
head(kk_readable)
```

### 4. Export Results

```r
write.csv(as.data.frame(kk_readable), 'kegg_results.csv', row.names = FALSE)
```

## KEGG vs GO

| Feature | KEGG | GO |
|---------|------|----|
| Focus | Pathways, reactions | Functions, processes |
| Structure | Curated pathway maps | Hierarchical ontology |
| Coverage | 4000+ organisms | Requires OrgDb |
| Update | Online database | OrgDb version |
| readable param | No (use setReadable) | Yes |

## Finding Organism Codes

```r
# Search for organism
search_kegg_organism('mouse')
#   kegg_code scientific_name    common_name
#        mmu   Mus musculus          mouse

search_kegg_organism('zebrafish')
#   kegg_code scientific_name    common_name
#        dre   Danio rerio      zebrafish
```

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

## KEGG Modules

KEGG modules are smaller functional units:

```r
mkk <- enrichMKEGG(
    gene = gene_list,
    organism = 'hsa'
)
```

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

## Visualization

See enrichment-visualization skill:

```r
# Dot plot
dotplot(kk, showCategory = 20)

# Bar plot
barplot(kk)

# Open pathway in browser
browseKEGG(kk, 'hsa04110')
```

## Common Issues

**No results found:**
- Check organism code is correct
- Verify gene IDs are Entrez format
- Increase pvalueCutoff

**Network error:**
- KEGG queries require internet
- Try `use_internal_data = TRUE` for cached data

**ID conversion fails:**
- Not all gene symbols map to Entrez
- Check conversion results for NAs
