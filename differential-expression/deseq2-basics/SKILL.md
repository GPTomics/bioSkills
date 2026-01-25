---
name: bio-de-deseq2-basics
description: Perform differential expression analysis using DESeq2 in R/Bioconductor. Use for analyzing RNA-seq count data, creating DESeqDataSet objects, running the DESeq workflow, and extracting results with log fold change shrinkage. Use when performing DE analysis with DESeq2.
tool_type: r
primary_tool: DESeq2
---

# DESeq2 Basics

Differential expression analysis using DESeq2 for RNA-seq count data.

## Required Libraries

```r
library(DESeq2)
library(apeglm)  # For lfcShrink with type='apeglm'
```

## Installation

```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')
BiocManager::install('DESeq2')
BiocManager::install('apeglm')
```

## Creating DESeqDataSet

### From Count Matrix

```r
# counts: matrix with genes as rows, samples as columns
# coldata: data frame with sample metadata (rownames must match colnames of counts)
dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData = coldata,
                               design = ~ condition)
```

### From SummarizedExperiment

```r
library(SummarizedExperiment)
dds <- DESeqDataSet(se, design = ~ condition)
```

### From tximport (Salmon/Kallisto)

```r
library(tximport)
txi <- tximport(files, type = 'salmon', tx2gene = tx2gene)
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)
```

## Standard DESeq2 Workflow

```r
# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData = coldata,
                               design = ~ condition)

# Pre-filter low count genes (recommended)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Set reference level for condition
dds$condition <- relevel(dds$condition, ref = 'control')

# Run DESeq2 pipeline (estimateSizeFactors, estimateDispersions, nbinomWaldTest)
dds <- DESeq(dds)

# Get results
res <- results(dds)

# Apply log fold change shrinkage (recommended for visualization/ranking)
resLFC <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')
```

## Design Formulas

```r
# Simple two-group comparison
design = ~ condition

# Controlling for batch effects
design = ~ batch + condition

# Interaction model
design = ~ genotype + treatment + genotype:treatment

# Multi-factor without interaction
design = ~ genotype + treatment
```

## Specifying Contrasts

```r
# See available coefficients
resultsNames(dds)

# Results by coefficient name
res <- results(dds, name = 'condition_treated_vs_control')

# Results by contrast (compare specific levels)
res <- results(dds, contrast = c('condition', 'treated', 'control'))

# Contrast with list format (for complex designs)
res <- results(dds, contrast = list('conditionB', 'conditionA'))
```

## Log Fold Change Shrinkage

```r
# apeglm method (default, recommended)
resLFC <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')

# ashr method (alternative)
resLFC <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'ashr')

# normal method (original, less recommended)
resLFC <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'normal')
```

## Setting Significance Thresholds

```r
# Default: padj < 0.1
res <- results(dds)

# Custom alpha threshold
res <- results(dds, alpha = 0.05)

# With log fold change threshold
res <- results(dds, lfcThreshold = 1)  # |log2FC| > 1
```

## Accessing DESeq2 Results

```r
# Summary of results
summary(res)

# Get significant genes
sig <- subset(res, padj < 0.05)

# Order by adjusted p-value
resOrdered <- res[order(res$padj),]

# Order by log fold change
resOrdered <- res[order(abs(res$log2FoldChange), decreasing = TRUE),]

# Convert to data frame
res_df <- as.data.frame(res)
```

## Result Columns

| Column | Description |
|--------|-------------|
| `baseMean` | Mean of normalized counts across all samples |
| `log2FoldChange` | Log2 fold change (treatment vs control) |
| `lfcSE` | Standard error of log2 fold change |
| `stat` | Wald statistic |
| `pvalue` | Raw p-value |
| `padj` | Adjusted p-value (Benjamini-Hochberg) |

## Normalization and Counts

```r
# Get normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# Get size factors
sizeFactors(dds)

# Variance stabilizing transformation (for visualization)
vsd <- vst(dds, blind = FALSE)

# Regularized log transformation (alternative, slower)
rld <- rlog(dds, blind = FALSE)
```

## Multi-Factor Designs

```r
# Design with batch correction
dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData = coldata,
                               design = ~ batch + condition)
dds <- DESeq(dds)

# Extract condition effect (controlling for batch)
res <- results(dds, name = 'condition_treated_vs_control')
```

## Interaction Models

```r
# Interaction between genotype and treatment
dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData = coldata,
                               design = ~ genotype + treatment + genotype:treatment)
dds <- DESeq(dds)

# Test interaction term
res_interaction <- results(dds, name = 'genotypeKO.treatmentdrug')

# Or use contrast for difference of differences
res_interaction <- results(dds, contrast = list(
    c('genotypeKO.treatmentdrug'),
    c()
))
```

## Likelihood Ratio Test

```r
# Compare full vs reduced model
dds <- DESeq(dds, test = 'LRT', reduced = ~ batch)

# Results from LRT
res <- results(dds)
```

## Pre-Filtering Strategies

```r
# Remove genes with low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Keep genes with at least n counts in at least k samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

# Filter by expression level
keep <- rowMeans(counts(dds, normalized = TRUE)) >= 10
dds <- dds[keep,]
```

## Working with Existing Objects

```r
# Update design formula
design(dds) <- ~ batch + condition
dds <- DESeq(dds)

# Subset samples
dds_subset <- dds[, dds$group == 'A']

# Subset genes
dds_genes <- dds[rownames(dds) %in% gene_list,]
```

## Exporting Results

```r
# Write to CSV
write.csv(as.data.frame(resOrdered), file = 'deseq2_results.csv')

# Write normalized counts
write.csv(as.data.frame(normalized_counts), file = 'normalized_counts.csv')
```

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "design matrix not full rank" | Confounded variables or missing levels | Check coldata for confounding |
| "counts matrix should be integers" | Non-integer counts (e.g., from tximport) | Use DESeqDataSetFromTximport() |
| "all samples have 0 counts" | Gene filtering issue | Check count matrix format |
| "factor levels not in colData" | Typo in design formula | Verify column names in coldata |

## Deprecated Features

| Feature | Status | Alternative |
|---------|--------|-------------|
| No-replicate designs | Removed (v1.22) | Require biological replicates |
| `betaPrior = TRUE` | Deprecated | Use `lfcShrink()` instead |
| `rlog()` for large datasets | Not recommended | Use `vst()` for >100 samples |

## Quick Reference: Workflow Steps

```r
# 1. Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ condition)

# 2. Pre-filter
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# 3. Set reference level
dds$condition <- relevel(dds$condition, ref = 'control')

# 4. Run DESeq2
dds <- DESeq(dds)

# 5. Get results with shrinkage
res <- lfcShrink(dds, coef = resultsNames(dds)[2], type = 'apeglm')

# 6. Filter significant genes
sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
```

## Related Skills

- edger-basics - Alternative DE analysis with edgeR
- de-visualization - MA plots, volcano plots, heatmaps
- de-results - Extract and export significant genes
