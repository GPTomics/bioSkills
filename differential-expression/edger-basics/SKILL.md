---
name: bio-de-edger-basics
description: Perform differential expression analysis using edgeR in R/Bioconductor. Use for analyzing RNA-seq count data with the quasi-likelihood F-test framework, creating DGEList objects, normalization, dispersion estimation, and statistical testing. Use when performing DE analysis with edgeR.
tool_type: r
primary_tool: edgeR
---

## Version Compatibility

Reference examples tested with: DESeq2 1.42+, edgeR 4.0+, limma 3.58+, scanpy 1.10+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

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

**Goal:** Construct an edgeR container from a count matrix with sample group information.

**Approach:** Wrap raw counts and group labels into a DGEList object for normalization and testing.

**"Load my RNA-seq counts into edgeR"** → Create a DGEList from a count matrix with sample group assignments and optional gene annotations.

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

**Goal:** Run the complete edgeR QL pipeline from raw counts to differentially expressed gene lists.

**Approach:** Filter, normalize (TMM), estimate dispersions, fit quasi-likelihood GLM, and test coefficients with the QL F-test.

**"Find differentially expressed genes between my groups"** → Test for significant expression differences using negative binomial models with quasi-likelihood F-tests.

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

**Goal:** Remove genes with insufficient expression to reduce noise and multiple testing burden.

**Approach:** Apply automatic or manual CPM/count thresholds requiring expression in a minimum number of samples.

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

**Goal:** Correct for differences in library composition between samples.

**Approach:** Compute TMM (or alternative) normalization factors that adjust effective library sizes.

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

**Goal:** Define the linear model structure for the experimental design.

**Approach:** Build model matrices encoding group, batch, and interaction terms for the GLM.

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

**Goal:** Estimate biological variability (dispersion) to parameterize the negative binomial model.

**Approach:** Compute common, trended, and gene-wise dispersions using empirical Bayes moderation.

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

**Goal:** Test for differential expression using the quasi-likelihood framework for robust inference.

**Approach:** Fit a QL GLM and test individual coefficients, contrasts, or multiple coefficients simultaneously.

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

**Goal:** Define specific pairwise or complex group comparisons for testing.

**Approach:** Use makeContrasts with a no-intercept design to specify arbitrary between-group differences.

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

**Goal:** Retrieve and filter DE results from the fitted model.

**Approach:** Use topTags to extract ranked gene lists with FDR-corrected p-values.

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

**Goal:** Perform a simple two-group DE test without a design matrix.

**Approach:** Use the classic edgeR exact test based on the negative binomial distribution.

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

**Goal:** Test for DE using likelihood ratio tests as an alternative to the QL F-test.

**Approach:** Fit a standard GLM and compare nested models via the likelihood ratio statistic.

```r
# Fit GLM
fit <- glmFit(y, design)

# Likelihood ratio test
lrt <- glmLRT(fit, coef = 2)
topTags(lrt)
```

## Treat Test (Log Fold Change Threshold)

**Goal:** Test whether genes exceed a minimum fold change threshold, not just differ from zero.

**Approach:** Use glmTreat to apply a fold change threshold directly in the statistical test.

```r
# Test for |logFC| > threshold
tr <- glmTreat(fit, coef = 2, lfc = log2(1.5))  # |FC| > 1.5
topTags(tr)
```

## Multi-Factor Designs

**Goal:** Test for condition effects while adjusting for batch or other covariates.

**Approach:** Include nuisance variables in the design matrix so the QL test controls for them.

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

**Goal:** Obtain normalized expression values for visualization and downstream analysis.

**Approach:** Compute CPM or log-CPM values using TMM-adjusted library sizes.

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

**Goal:** Save DE results and normalized counts to files for sharing or downstream tools.

**Approach:** Extract all results via topTags and write as CSV alongside CPM values.

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

| Scenario | Recommended | Rationale |
|----------|-------------|-----------|
| Standard bulk RNA-seq (n >= 3/group) | Either — expect ~90% overlap | Both perform well; concordance is high |
| Very small sample size (n = 2-3/group) | edgeR QL F-test | QL framework provides tighter FPR control with few replicates |
| Large datasets (>100 samples, many conditions) | edgeR | Faster; C++ backend in v4 |
| Salmon/Kallisto quantification | DESeq2 via tximport | DESeqDataSetFromTximport handles offset matrix natively |
| Need LFC shrinkage for ranking/GSEA | DESeq2 | apeglm/ashr shrinkage built in; edgeR has no equivalent |
| Formal fold-change threshold testing | edgeR `glmTreat` | More flexible than DESeq2 `lfcThreshold` |
| Python-only environment | DESeq2 (via PyDESeq2) | No edgeR Python equivalent |
| Omnibus test (>= 3 groups) | Either (LRT) | Both support likelihood ratio tests |

If overlap between DESeq2 and edgeR is < 60%, investigate differences in filtering, normalization, or dispersion estimation — discordance usually indicates a modeling issue, not a tool difference.

## Decision Guidance

### Test Selection

| Test | When to Use | Key Property |
|------|-------------|--------------|
| QL F-test (`glmQLFit` + `glmQLFTest`) | **Default for most experiments** | Best FPR control with small n; accounts for dispersion uncertainty |
| LRT (`glmFit` + `glmLRT`) | Large samples (n >= 6 per group); complex designs where QL is too conservative | More powerful but can be anti-conservative with few replicates |
| Exact test (`exactTest`) | Simple two-group comparison only | No design matrix; cannot adjust for covariates |
| TREAT (`glmTreat`) | Testing against a minimum fold-change threshold | Tests H0: \|LFC\| <= threshold, not H0: LFC = 0 |

QL F-test p-values are always >= LRT p-values. In null comparisons (replicates vs replicates), QL consistently returns ~0 false positives while LRT can return many.

### edgeR v4 Changes

- Constant NB dispersion estimated from the most highly expressed genes (v3 used trended dispersions)
- `estimateDisp()` is now optional before `glmQLFit()` (but still needed for BCV plots)
- Fractional count support for transcript quantification uncertainty (Gibbs sampling)
- C++ backend for model fitting (faster)
- `decidetestsDGE()` removed; use `decideTests()` instead
- `mglmLS()`, `mglmSimple()` retired; use `mglmLevenberg()` or `mglmOneWay()`

### filterByExpr Internals

Default parameters: `min.count = 10`, `min.total.count = 15`, `large.n = 10`, `min.prop = 0.7`.

Algorithm:
1. Convert `min.count` to CPM cutoff: `min.count / median(lib.size) * 1e6`
2. Determine minimum number of samples `n` from design (smallest group size)
3. If group size > `large.n`: `n = large.n + (group_size - large.n) * min.prop`
4. Keep gene if CPM >= cutoff in >= n samples AND total count >= `min.total.count`

Example: median library 51M reads -> CPM cutoff ~0.2; 3 replicates per group -> need CPM >= 0.2 in >= 3 samples.

### Normalization Caveats

TMM/RLE both assume most genes are NOT DE. This assumption breaks in:
- Prokaryotic stress responses (majority of genes may be DE)
- Extreme perturbations or cross-species comparisons
- Single-cell with many zeros

When violated, consider spike-in normalization or supplying known stable reference genes. For prokaryotic RNA-seq, also note: use non-spliced aligners (BWA-MEM, Bowtie2), rRNA depletion is essential, and KEGG organism codes are strain-specific.

## Related Skills

- deseq2-basics - Alternative DE analysis with DESeq2
- de-visualization - MA plots, volcano plots, heatmaps
- de-results - Extract and export significant genes
