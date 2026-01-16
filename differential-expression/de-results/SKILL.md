---
name: bio-de-results
description: Extract, filter, annotate, and export differential expression results from DESeq2 or edgeR. Use for identifying significant genes, applying multiple testing corrections, adding gene annotations, and preparing results for downstream analysis.
tool_type: r
primary_tool: DESeq2
---

# DE Results

Extract, filter, and export differential expression results.

## Required Libraries

```r
library(DESeq2)  # or library(edgeR)
library(dplyr)   # For data manipulation
```

## Extracting DESeq2 Results

```r
# Basic results
res <- results(dds)

# With specific alpha (adjusted p-value threshold)
res <- results(dds, alpha = 0.05)

# With log fold change shrinkage
res <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')

# Convert to data frame
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
```

## Extracting edgeR Results

```r
# Get all results
results <- topTags(qlf, n = Inf)$table

# Add gene column
results$gene <- rownames(results)
```

## Filtering Significant Genes

### By Adjusted P-value

```r
# DESeq2
sig_genes <- subset(res, padj < 0.05)

# edgeR
sig_genes <- subset(results, FDR < 0.05)

# Using dplyr
sig_genes <- res_df %>%
    filter(padj < 0.05) %>%
    arrange(padj)
```

### By Fold Change

```r
# Absolute log2 fold change > 1 (2-fold change)
sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)

# Up-regulated only
up_genes <- subset(res, padj < 0.05 & log2FoldChange > 1)

# Down-regulated only
down_genes <- subset(res, padj < 0.05 & log2FoldChange < -1)
```

### Combined Filters

```r
# Stringent filtering
sig_genes <- res_df %>%
    filter(padj < 0.01,
           abs(log2FoldChange) > 1,
           baseMean > 10) %>%
    arrange(padj)
```

## Ordering Results

```r
# By adjusted p-value (most significant first)
res_ordered <- res[order(res$padj), ]

# By absolute fold change (largest changes first)
res_ordered <- res[order(abs(res$log2FoldChange), decreasing = TRUE), ]

# By base mean expression
res_ordered <- res[order(res$baseMean, decreasing = TRUE), ]

# Combined: significant genes ordered by fold change
sig_ordered <- res_df %>%
    filter(padj < 0.05) %>%
    arrange(desc(abs(log2FoldChange)))
```

## Summary Statistics

```r
# DESeq2 summary
summary(res)

# Manual counts
n_tested <- sum(!is.na(res$padj))
n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
n_up <- sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
n_down <- sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)

cat(sprintf('Tested: %d genes\n', n_tested))
cat(sprintf('Significant (padj < 0.05): %d genes\n', n_sig))
cat(sprintf('Up-regulated: %d genes\n', n_up))
cat(sprintf('Down-regulated: %d genes\n', n_down))

# edgeR summary
summary(decideTests(qlf))
```

## Adding Gene Annotations

### From Bioconductor Annotation Package

```r
library(org.Hs.eg.db)  # Human; use org.Mm.eg.db for mouse

# If gene IDs are Ensembl
res_df$symbol <- mapIds(org.Hs.eg.db,
                         keys = rownames(res_df),
                         column = 'SYMBOL',
                         keytype = 'ENSEMBL',
                         multiVals = 'first')

res_df$entrez <- mapIds(org.Hs.eg.db,
                         keys = rownames(res_df),
                         column = 'ENTREZID',
                         keytype = 'ENSEMBL',
                         multiVals = 'first')

res_df$description <- mapIds(org.Hs.eg.db,
                              keys = rownames(res_df),
                              column = 'GENENAME',
                              keytype = 'ENSEMBL',
                              multiVals = 'first')
```

### From BioMart

```r
library(biomaRt)

mart <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')

annotations <- getBM(
    attributes = c('ensembl_gene_id', 'external_gene_name', 'description'),
    filters = 'ensembl_gene_id',
    values = rownames(res_df),
    mart = mart
)

# Merge with results
res_annotated <- merge(res_df, annotations,
                        by.x = 'row.names', by.y = 'ensembl_gene_id',
                        all.x = TRUE)
```

### From Custom File

```r
# Load annotation file
gene_info <- read.csv('gene_annotations.csv')

# Merge with results
res_annotated <- merge(res_df, gene_info, by = 'gene', all.x = TRUE)
```

## Exporting Results

### To CSV

```r
# All results
write.csv(res_df, file = 'deseq2_all_results.csv', row.names = FALSE)

# Significant only
sig_genes <- res_df %>% filter(padj < 0.05)
write.csv(sig_genes, file = 'deseq2_significant.csv', row.names = FALSE)
```

### To Excel

```r
library(openxlsx)

# Create workbook with multiple sheets
wb <- createWorkbook()

addWorksheet(wb, 'All Results')
writeData(wb, 'All Results', res_df)

addWorksheet(wb, 'Significant')
writeData(wb, 'Significant', sig_genes)

addWorksheet(wb, 'Up-regulated')
writeData(wb, 'Up-regulated', up_genes)

addWorksheet(wb, 'Down-regulated')
writeData(wb, 'Down-regulated', down_genes)

saveWorkbook(wb, 'de_results.xlsx', overwrite = TRUE)
```

### Gene Lists for Pathway Analysis

```r
# Just gene IDs for GO/KEGG analysis
sig_gene_list <- rownames(subset(res, padj < 0.05))
write.table(sig_gene_list, file = 'significant_genes.txt',
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# With fold changes for GSEA
gsea_input <- res_df %>%
    filter(!is.na(log2FoldChange)) %>%
    select(gene, log2FoldChange) %>%
    arrange(desc(log2FoldChange))
write.table(gsea_input, file = 'gsea_input.rnk',
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
```

## Comparing Results Between Methods

```r
# Get significant genes from both methods
deseq2_sig <- rownames(subset(deseq2_res, padj < 0.05))
edger_sig <- rownames(subset(edger_results, FDR < 0.05))

# Overlap
common <- intersect(deseq2_sig, edger_sig)
deseq2_only <- setdiff(deseq2_sig, edger_sig)
edger_only <- setdiff(edger_sig, deseq2_sig)

cat(sprintf('DESeq2 significant: %d\n', length(deseq2_sig)))
cat(sprintf('edgeR significant: %d\n', length(edger_sig)))
cat(sprintf('Common: %d\n', length(common)))
cat(sprintf('DESeq2 only: %d\n', length(deseq2_only)))
cat(sprintf('edgeR only: %d\n', length(edger_only)))

# Venn diagram
library(VennDiagram)
venn.diagram(
    x = list(DESeq2 = deseq2_sig, edgeR = edger_sig),
    filename = 'de_overlap.png',
    fill = c('steelblue', 'coral')
)
```

## Multiple Testing Correction

```r
# DESeq2 uses Benjamini-Hochberg by default
# To use different methods:

# Independent Hypothesis Weighting (more powerful)
library(IHW)
res_ihw <- results(dds, filterFun = ihw)

# Manual p-value adjustment
res_df$padj_bonferroni <- p.adjust(res_df$pvalue, method = 'bonferroni')
res_df$padj_bh <- p.adjust(res_df$pvalue, method = 'BH')
res_df$padj_fdr <- p.adjust(res_df$pvalue, method = 'fdr')
```

## Handling NA Values

```r
# Count NAs
sum(is.na(res$padj))

# Remove genes with NA padj
res_complete <- res[!is.na(res$padj), ]

# Understand why NAs occur
# - baseMean = 0: No counts
# - NA only in padj: Outlier or low count filtered by independent filtering

# Check outliers
res[which(is.na(res$pvalue) & res$baseMean > 0), ]
```

## Quick Reference: Result Columns

### DESeq2

| Column | Description |
|--------|-------------|
| `baseMean` | Mean normalized counts |
| `log2FoldChange` | Log2 fold change |
| `lfcSE` | Standard error of LFC |
| `stat` | Wald statistic |
| `pvalue` | Raw p-value |
| `padj` | Adjusted p-value (BH) |

### edgeR

| Column | Description |
|--------|-------------|
| `logFC` | Log2 fold change |
| `logCPM` | Average log2 CPM |
| `F` | Quasi-likelihood F-statistic |
| `PValue` | Raw p-value |
| `FDR` | False discovery rate |

## Related Skills

- **deseq2-basics** - Run DESeq2 analysis
- **edger-basics** - Run edgeR analysis
- **de-visualization** - Visualize results
- **pathway-analysis** - Use gene lists for enrichment
