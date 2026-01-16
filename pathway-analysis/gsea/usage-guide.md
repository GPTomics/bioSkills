# GSEA Usage Guide

Gene Set Enrichment Analysis (GSEA) tests whether genes in predefined sets show coordinated changes across conditions.

## GSEA vs Over-Representation

| Feature | Over-Representation | GSEA |
|---------|---------------------|------|
| Input | Gene list (significant only) | Ranked gene list (all genes) |
| Cutoff | Requires significance threshold | No arbitrary cutoff |
| Detection | Strong individual changes | Coordinated subtle changes |
| Functions | enrichGO, enrichKEGG | gseGO, gseKEGG |

## When to Use GSEA

- You want to use information from all genes
- You're looking for coordinated, subtle changes
- Over-representation finds few or no enriched terms
- You have good expression statistics for all genes

## Prerequisites

```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install(c('clusterProfiler', 'org.Hs.eg.db'))

# For MSigDB gene sets:
install.packages('msigdbr')
```

## Basic Workflow

### 1. Create Ranked Gene List

The most critical step - genes must be:
- Named (gene IDs)
- Numeric (statistic)
- Sorted (decreasing order)

```r
library(clusterProfiler)
library(org.Hs.eg.db)

de_results <- read.csv('deseq2_results.csv')

# Use log2FoldChange as ranking statistic
gene_list <- de_results$log2FoldChange
names(gene_list) <- de_results$gene_id

# Remove NAs and sort
gene_list <- gene_list[!is.na(gene_list)]
gene_list <- sort(gene_list, decreasing = TRUE)
```

### 2. Convert to Entrez IDs

```r
# Convert symbols to Entrez
gene_ids <- bitr(names(gene_list), fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)

# Create Entrez-named list
gene_list_entrez <- gene_list[names(gene_list) %in% gene_ids$SYMBOL]
names(gene_list_entrez) <- gene_ids$ENTREZID[match(names(gene_list_entrez), gene_ids$SYMBOL)]
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)
```

### 3. Run GSEA

```r
# GO GSEA
gse <- gseGO(
    geneList = gene_list_entrez,
    OrgDb = org.Hs.eg.db,
    ont = 'BP',
    pvalueCutoff = 0.05,
    verbose = FALSE
)
```

### 4. View Results

```r
head(gse)
results <- as.data.frame(gse)
```

## Choosing a Ranking Statistic

| Statistic | Formula | Best For |
|-----------|---------|----------|
| log2FC | log2FoldChange | Magnitude of change |
| Signed p | -log10(p) * sign(FC) | Both significance and direction |
| Wald | stat column from DESeq2 | Pre-computed statistic |
| t-statistic | From limma | Moderated statistics |

Signed p-value often works best:
```r
gene_list <- -log10(de_results$pvalue) * sign(de_results$log2FoldChange)
```

## Understanding Results

| Column | Description |
|--------|-------------|
| ID | Gene set ID |
| Description | Gene set name |
| setSize | Genes in set |
| enrichmentScore | Raw ES |
| NES | Normalized Enrichment Score |
| pvalue | Nominal p-value |
| p.adjust | FDR-adjusted p-value |
| core_enrichment | Leading edge genes |

### Interpreting NES

- **Positive NES**: Gene set genes tend to be upregulated
- **Negative NES**: Gene set genes tend to be downregulated
- **|NES| > 1.5**: Strong enrichment
- **|NES| > 2.0**: Very strong enrichment

## Visualization

See enrichment-visualization skill:

```r
# GSEA plot (running score)
gseaplot2(gse, geneSetID = 1, title = gse$Description[1])

# Ridge plot (multiple gene sets)
ridgeplot(gse)

# Dot plot
dotplot(gse, showCategory = 20)
```

## Using MSigDB Gene Sets

```r
library(msigdbr)

# Get Hallmark gene sets
hallmarks <- msigdbr(species = 'Homo sapiens', category = 'H')
hallmarks_t2g <- hallmarks[, c('gs_name', 'entrez_gene')]

# Run GSEA
gse <- GSEA(geneList = gene_list_entrez, TERM2GENE = hallmarks_t2g)
```

## Common Issues

**Error: geneList must be sorted:**
- Ensure `sort(gene_list, decreasing = TRUE)`

**No enriched terms:**
- Try different ranking statistic
- Increase pvalueCutoff
- Check gene ID conversion

**Many NAs in ranked list:**
- Remove NAs before sorting
- Check for missing values in DE results
