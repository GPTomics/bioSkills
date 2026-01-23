---
name: bio-pathway-kegg-pathways
description: KEGG pathway and module enrichment analysis using clusterProfiler enrichKEGG and enrichMKEGG. Use when identifying metabolic and signaling pathways over-represented in a gene list. Supports 4000+ organisms via KEGG online database.
tool_type: r
primary_tool: clusterProfiler
---

# KEGG Pathway Enrichment

## Core Pattern

```r
library(clusterProfiler)

kk <- enrichKEGG(
    gene = gene_list,           # Character vector of gene IDs
    organism = 'hsa',           # KEGG organism code
    pvalueCutoff = 0.05,
    pAdjustMethod = 'BH'
)
```

## Prepare Gene List

```r
library(org.Hs.eg.db)

de_results <- read.csv('de_results.csv')
sig_genes <- de_results$gene_id[de_results$padj < 0.05 & abs(de_results$log2FoldChange) > 1]

# KEGG requires NCBI Entrez gene IDs (kegg, ncbi-geneid)
gene_ids <- bitr(sig_genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
gene_list <- gene_ids$ENTREZID
```

## KEGG ID Conversion

```r
# Convert between KEGG and other IDs
kegg_ids <- bitr_kegg(gene_list, fromType = 'ncbi-geneid', toType = 'kegg', organism = 'hsa')

# Available types: kegg, ncbi-geneid, ncbi-proteinid, uniprot
```

## Run KEGG Pathway Enrichment

```r
kk <- enrichKEGG(
    gene = gene_list,
    organism = 'hsa',
    keyType = 'ncbi-geneid',    # or 'kegg'
    pvalueCutoff = 0.05,
    pAdjustMethod = 'BH',
    minGSSize = 10,
    maxGSSize = 500
)

# View results
head(kk)
results <- as.data.frame(kk)
```

## Make Results Readable

```r
# enrichKEGG does NOT have readable parameter - use setReadable
library(org.Hs.eg.db)
kk_readable <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
```

## KEGG Module Enrichment

```r
# KEGG modules are smaller functional units than pathways
mkk <- enrichMKEGG(
    gene = gene_list,
    organism = 'hsa',
    pvalueCutoff = 0.05
)
```

## Common Organism Codes

| Organism | Code | Common Name |
|----------|------|-------------|
| hsa | Human | Homo sapiens |
| mmu | Mouse | Mus musculus |
| rno | Rat | Rattus norvegicus |
| dre | Zebrafish | Danio rerio |
| dme | Fruit fly | Drosophila melanogaster |
| cel | Worm | C. elegans |
| sce | Yeast | S. cerevisiae |
| ath | Arabidopsis | A. thaliana |
| eco | E. coli K-12 | |

```r
# Find organism codes
search_kegg_organism('mouse')
search_kegg_organism('zebrafish')
```

## With Background Universe

```r
all_genes <- de_results$gene_id
universe_ids <- bitr(all_genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)

kk <- enrichKEGG(
    gene = gene_list,
    universe = universe_ids$ENTREZID,
    organism = 'hsa',
    pvalueCutoff = 0.05
)
```

## Extract and Export Results

```r
# Convert to data frame
results_df <- as.data.frame(kk)

# Key columns: ID (pathway), Description, GeneRatio, BgRatio, pvalue, p.adjust, geneID, Count

# Export
write.csv(results_df, 'kegg_enrichment_results.csv', row.names = FALSE)

# Get genes in a specific pathway
pathway_genes <- kk@geneSets[['hsa04110']]  # Cell cycle
```

## Browse KEGG Pathways

```r
# View pathway in browser (opens KEGG website)
browseKEGG(kk, 'hsa04110')

# Download pathway image
library(pathview)
pathview(gene.data = gene_list, pathway.id = 'hsa04110', species = 'hsa')
```

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| gene | required | Vector of gene IDs |
| organism | hsa | KEGG organism code |
| keyType | kegg | Input ID type |
| pvalueCutoff | 0.05 | P-value threshold |
| qvalueCutoff | 0.2 | Q-value threshold |
| pAdjustMethod | BH | Adjustment method |
| universe | NULL | Background genes |
| minGSSize | 10 | Min genes per pathway |
| maxGSSize | 500 | Max genes per pathway |
| use_internal_data | FALSE | Use local KEGG data |

## Compare Multiple Gene Lists

```r
# Compare KEGG enrichment across groups
gene_lists <- list(
    up = up_genes,
    down = down_genes
)

ck <- compareCluster(
    geneClusters = gene_lists,
    fun = 'enrichKEGG',
    organism = 'hsa'
)

dotplot(ck)
```

## Notes

- **No readable parameter** - use `setReadable()` with OrgDb
- **Requires internet** - queries KEGG database online
- **use_internal_data** - set TRUE to use cached KEGG data (may be outdated)
- **Pathway IDs** - format is organism code + 5 digits (e.g., hsa04110)

## Related Skills

- go-enrichment - Gene Ontology enrichment analysis
- gsea - GSEA using KEGG pathways (gseKEGG)
- enrichment-visualization - Visualize KEGG results
