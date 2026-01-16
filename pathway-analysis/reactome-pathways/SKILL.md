---
name: bio-pathway-reactome
description: Reactome pathway enrichment using ReactomePA package. Performs over-representation analysis and GSEA on Reactome pathways, a curated peer-reviewed knowledgebase of biological pathways. Supports visualization and pathway hierarchy exploration.
tool_type: r
primary_tool: ReactomePA
---

# Reactome Pathway Enrichment

## Core Pattern - Over-Representation Analysis

```r
library(ReactomePA)
library(org.Hs.eg.db)

pathway_result <- enrichPathway(
    gene = entrez_ids,         # Character vector of Entrez IDs
    organism = 'human',        # human, rat, mouse, celegans, yeast, zebrafish, fly
    pvalueCutoff = 0.05,
    pAdjustMethod = 'BH',
    readable = TRUE            # Convert to gene symbols
)

head(as.data.frame(pathway_result))
```

## Prepare Gene List from DE Results

```r
library(clusterProfiler)

de_results <- read.csv('de_results.csv')
sig_genes <- de_results[de_results$padj < 0.05 & abs(de_results$log2FoldChange) > 1, 'gene_symbol']

gene_ids <- bitr(sig_genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
entrez_ids <- gene_ids$ENTREZID
```

## GSEA on Reactome Pathways

```r
# Create ranked gene list (named vector sorted by statistic)
gene_list <- de_results$log2FoldChange
names(gene_list) <- de_results$entrez_id
gene_list <- sort(gene_list, decreasing = TRUE)

gsea_result <- gsePathway(
    geneList = gene_list,
    organism = 'human',
    pvalueCutoff = 0.05,
    pAdjustMethod = 'BH',
    verbose = FALSE
)

head(as.data.frame(gsea_result))
```

## With Background Universe

```r
all_genes <- de_results$entrez_id  # All tested genes

pathway_result <- enrichPathway(
    gene = entrez_ids,
    universe = all_genes,      # Background gene set
    organism = 'human',
    pvalueCutoff = 0.05,
    readable = TRUE
)
```

## Visualization

```r
library(enrichplot)

# Dot plot
dotplot(pathway_result, showCategory = 15)

# Bar plot
barplot(pathway_result, showCategory = 15)

# Enrichment map (requires pairwise_termsim first)
pathway_result <- pairwise_termsim(pathway_result)
emapplot(pathway_result)

# Gene-concept network
cnetplot(pathway_result, categorySize = 'pvalue')

# GSEA plot
gseaplot2(gsea_result, geneSetID = 1:3)
```

## View Pathway in Browser

```r
# Open pathway in Reactome browser
viewPathway('R-HSA-109582', organism = 'human')  # Uses pathway ID

# Get pathway ID from results
top_pathway_id <- pathway_result@result$ID[1]
viewPathway(top_pathway_id, organism = 'human')
```

## Export Results

```r
results_df <- as.data.frame(pathway_result)
write.csv(results_df, 'reactome_enrichment.csv', row.names = FALSE)

# Key columns: ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, geneID, Count
```

## Different Organisms

```r
# Mouse
pathway_mouse <- enrichPathway(gene = mouse_entrez, organism = 'mouse', readable = TRUE)

# Rat
pathway_rat <- enrichPathway(gene = rat_entrez, organism = 'rat', readable = TRUE)

# Zebrafish
pathway_zfish <- enrichPathway(gene = zfish_entrez, organism = 'zebrafish', readable = TRUE)

# Supported: human, rat, mouse, celegans, yeast, zebrafish, fly
```

## Compare Clusters

```r
# Compare pathways across multiple gene lists
gene_clusters <- list(
    upregulated = up_genes,
    downregulated = down_genes
)

compare_result <- compareCluster(
    geneClusters = gene_clusters,
    fun = 'enrichPathway',
    organism = 'human',
    pvalueCutoff = 0.05
)

dotplot(compare_result)
```

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| gene | required | Vector of Entrez IDs |
| organism | human | Species name |
| pvalueCutoff | 0.05 | P-value threshold |
| pAdjustMethod | BH | Adjustment method |
| universe | NULL | Background genes |
| minGSSize | 10 | Min genes per pathway |
| maxGSSize | 500 | Max genes per pathway |
| readable | FALSE | Convert to symbols |

## Supported Organisms

| Organism | Name | OrgDb |
|----------|------|-------|
| Human | human | org.Hs.eg.db |
| Mouse | mouse | org.Mm.eg.db |
| Rat | rat | org.Rn.eg.db |
| Zebrafish | zebrafish | org.Dr.eg.db |
| Fly | fly | org.Dm.eg.db |
| C. elegans | celegans | org.Ce.eg.db |
| Yeast | yeast | org.Sc.sgd.db |

## Related Skills

- go-enrichment - Gene Ontology enrichment
- kegg-pathways - KEGG pathway enrichment
- wikipathways - WikiPathways enrichment
- gsea - Gene Set Enrichment Analysis
- enrichment-visualization - Visualization functions
