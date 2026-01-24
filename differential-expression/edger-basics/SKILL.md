---
name: bio-de-edger-basics
description: Perform differential expression analysis using edgeR in R/Bioconductor. Use for analyzing RNA-seq count data with the quasi-likelihood F-test framework, creating DGEList objects, normalization, dispersion estimation, and statistical testing. Use when performing DE analysis with edgeR.
tool_type: r
primary_tool: edgeR
---

# edgeR Basics

Differential expression analysis using edgeR's quasi-likelihood framework for RNA-seq count data.

## Required Libraries

```r
library(edgeR)
library(limma)  # For design matrices and voom
```

## Installation

```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')
BiocManager::install('edgeR')
```

## Creating DGEList Object

```r
# From count matrix
# counts: matrix with genes as rows, samples as columns
# group: factor indicating sample groups
y <- DGEList(counts = counts, group = group)

# With gene annotation
y <- DGEList(counts = counts, group = group, genes = gene_info)

# Check structure
y
```

## Standard edgeR Workflow (Quasi-Likelihood)

```r
# Create DGEList
y <- DGEList(counts = counts, group = group)

# Filter low-expression genes
keep <- filterByExpr(y, group = group)
y <- y[keep, , keep.lib.sizes = FALSE]

# Normalize (TMM by default)
y <- calcNormFactors(y)

# Create design matrix
design <- model.matrix(~ group)

# Estimate dispersion (optional in edgeR v4+ but improves BCV plots)
y <- estimateDisp(y, design)

# Fit quasi-likelihood model
fit <- glmQLFit(y, design)

# Perform quasi-likelihood F-test
qlf <- glmQLFTest(fit, coef = 2)

# View top genes
topTags(qlf)
```

## Filtering Low-Expression Genes

```r
# Automatic filtering (recommended)
keep <- filterByExpr(y, group = group)
y <- y[keep, , keep.lib.sizes = FALSE]

# Manual filtering: CPM threshold
keep <- rowSums(cpm(y) > 1) >= 2  # At least 2 samples with CPM > 1
y <- y[keep, , keep.lib.sizes = FALSE]

# Filter by minimum counts
keep <- rowSums(y$counts >= 10) >= 3  # At least 3 samples with 10+ counts
y <- y[keep, , keep.lib.sizes = FALSE]
```

## Normalization Methods

```r
# TMM normalization (default, recommended)
y <- calcNormFactors(y, method = 'TMM')

# Alternative methods
y <- calcNormFactors(y, method = 'RLE')      # Relative Log Expression
y <- calcNormFactors(y, method = 'upperquartile')
y <- calcNormFactors(y, method = 'none')     # No normalization

# View normalization factors
y$samples$norm.factors
```

## Design Matrices

```r
# Simple two-group comparison
design <- model.matrix(~ group)

# With batch correction
design <- model.matrix(~ batch + group)

# Interaction model
design <- model.matrix(~ genotype + treatment + genotype:treatment)

# No intercept (for direct group comparisons)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
```

## Dispersion Estimation

```r
# Estimate all dispersions
y <- estimateDisp(y, design)

# Or estimate separately
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

# View dispersions
y$common.dispersion
y$trended.dispersion
y$tagwise.dispersion

# Plot BCV (biological coefficient of variation)
plotBCV(y)
```

## Quasi-Likelihood Testing

```r
# Fit QL model
fit <- glmQLFit(y, design)

# Test specific coefficient
qlf <- glmQLFTest(fit, coef = 2)

# Test with contrast
contrast <- makeContrasts(groupB - groupA, levels = design)
qlf <- glmQLFTest(fit, contrast = contrast)

# Test multiple coefficients (ANOVA-like)
qlf <- glmQLFTest(fit, coef = 2:3)
```

## Making Contrasts

```r
# Design without intercept
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

# Pairwise comparisons
contrast <- makeContrasts(
    TreatedVsControl = treated - control,
    DrugAVsControl = drugA - control,
    DrugBVsControl = drugB - control,
    DrugAVsDrugB = drugA - drugB,
    levels = design
)

# Test each contrast
qlf_treated <- glmQLFTest(fit, contrast = contrast[, 'TreatedVsControl'])
qlf_drugA <- glmQLFTest(fit, contrast = contrast[, 'DrugAVsControl'])
```

## Accessing Results

```r
# Top differentially expressed genes
topTags(qlf, n = 20)

# All results as data frame
results <- topTags(qlf, n = Inf)$table

# Summary of DE genes at different thresholds
summary(decideTests(qlf))

# Get DE genes with specific cutoffs
de_genes <- topTags(qlf, n = Inf, p.value = 0.05)$table
```

## Result Columns

| Column | Description |
|--------|-------------|
| `logFC` | Log2 fold change |
| `logCPM` | Average log2 counts per million |
| `F` | Quasi-likelihood F-statistic |
| `PValue` | Raw p-value |
| `FDR` | False discovery rate (adjusted p-value) |

## Alternative: Exact Test (Classic edgeR)

```r
# For simple two-group comparison only
y <- DGEList(counts = counts, group = group)
y <- calcNormFactors(y)
y <- estimateDisp(y)

# Exact test
et <- exactTest(y)
topTags(et)
```

## Alternative: glmLRT (Likelihood Ratio Test)

```r
# Fit GLM
fit <- glmFit(y, design)

# Likelihood ratio test
lrt <- glmLRT(fit, coef = 2)
topTags(lrt)
```

## Treat Test (Log Fold Change Threshold)

```r
# Test for |logFC| > threshold
tr <- glmTreat(fit, coef = 2, lfc = log2(1.5))  # |FC| > 1.5
topTags(tr)
```

## Multi-Factor Designs

```r
# Design with batch correction
design <- model.matrix(~ batch + condition, data = sample_info)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

# Test condition effect (controlling for batch)
# Condition coefficient is typically the last
qlf <- glmQLFTest(fit, coef = ncol(design))
```

## Getting Normalized Counts

```r
# Counts per million (CPM)
cpm_values <- cpm(y)

# Log2 CPM
log_cpm <- cpm(y, log = TRUE)

# RPM (reads per million, same as CPM)
rpm_values <- cpm(y)

# With prior count for log transformation
log_cpm <- cpm(y, log = TRUE, prior.count = 2)
```

## Exporting Results

```r
# Get all results
all_results <- topTags(qlf, n = Inf)$table

# Add gene IDs as column
all_results$gene_id <- rownames(all_results)

# Write to file
write.csv(all_results, file = 'edger_results.csv', row.names = FALSE)

# Export normalized counts
write.csv(cpm(y), file = 'cpm_values.csv')
```

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "design matrix not full rank" | Confounded variables | Check sample metadata |
| "No residual df" | Too few samples | Need more replicates |
| "NA/NaN/Inf" | Zero counts in all samples | Filter more stringently |

## Deprecated/Changed Functions

| Old | Status | New |
|-----|--------|-----|
| `decidetestsDGE()` | Removed (v4.4) | `decideTests()` |
| `glmFit()` + `glmLRT()` | Still works | Prefer `glmQLFit()` + `glmQLFTest()` |
| `estimateDisp()` | Optional (v4+) | `glmQLFit()` estimates internally |
| `mglmLS()`, `mglmSimple()` | Retired | `mglmLevenberg()` or `mglmOneWay()` |

**Note:** `calcNormFactors()` and `normLibSizes()` are synonyms - both work.

## Quick Reference: Workflow Steps

```r
# 1. Create DGEList
y <- DGEList(counts = counts, group = group)

# 2. Filter low-expression genes
keep <- filterByExpr(y, group = group)
y <- y[keep, , keep.lib.sizes = FALSE]

# 3. Normalize
y <- calcNormFactors(y)

# 4. Create design matrix
design <- model.matrix(~ group)

# 5. Estimate dispersion (optional in v4+)
y <- estimateDisp(y, design)

# 6. Fit quasi-likelihood model
fit <- glmQLFit(y, design)

# 7. Test for DE
qlf <- glmQLFTest(fit, coef = 2)

# 8. Get results
topTags(qlf, n = 20)
```

## Choosing edgeR vs DESeq2

| Aspect | edgeR | DESeq2 |
|--------|-------|--------|
| Model | Negative binomial + QL | Negative binomial |
| Shrinkage | Empirical Bayes on dispersions | Shrinkage estimators for LFC |
| Small samples | Robust with QL framework | Good with shrinkage |
| Speed | Generally faster | Slower for large datasets |
| Output | F-statistic, FDR | Wald statistic, padj |

## Related Skills

- deseq2-basics - Alternative DE analysis with DESeq2
- de-visualization - MA plots, volcano plots, heatmaps
- de-results - Extract and export significant genes
