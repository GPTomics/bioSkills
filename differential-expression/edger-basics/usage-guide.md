# edgeR Basics - Usage Guide

## Overview

This skill covers differential expression analysis using edgeR, a powerful Bioconductor package that uses empirical Bayes methods to moderate gene-wise dispersion estimates. The quasi-likelihood F-test framework provides robust results even with small sample sizes.

## Prerequisites

```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')
BiocManager::install('edgeR')
```

## Quick Start

Tell your AI agent what you want to do:
- "Run edgeR analysis on my count matrix with treated vs control comparison"
- "Perform differential expression with edgeR controlling for batch"
- "Get top differentially expressed genes with FDR < 0.05"

## Example Prompts

### Basic Analysis
> "Create a DGEList from my count matrix"

> "Run the standard edgeR quasi-likelihood workflow"

> "Normalize my RNA-seq data with TMM"

### Design and Contrasts
> "Set up edgeR with multiple comparisons"

> "Create contrasts to compare all treatments against control"

> "Test for interaction effects between genotype and treatment"

### Results
> "Get genes with FDR < 0.05 and |logFC| > 1"

> "Export edgeR results to CSV"

> "Compare edgeR and DESeq2 results"

## What the Agent Will Do

1. Create DGEList object from counts
2. Filter low-expression genes with filterByExpr()
3. Normalize using TMM (calcNormFactors)
4. Set up appropriate design matrix
5. Estimate dispersions (estimateDisp)
6. Fit quasi-likelihood model (glmQLFit)
7. Perform statistical testing (glmQLFTest)
8. Extract and filter results

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| Count matrix | Integer matrix | Genes (rows) x Samples (columns) |
| Group factor | Factor | Sample group assignments |
| Design matrix | Matrix | From model.matrix() |

## Tips

- Use `filterByExpr()` for automatic low-count filtering -- skipping it inflates the multiple-testing burden; edgeR has no automatic independent filtering like DESeq2
- The quasi-likelihood framework (glmQLFit) is recommended over exact test; QL F-test is the modern default with proper finite-sample FPR control
- Set BOTH `estimateDisp(robust=TRUE)` AND `glmQLFit(robust=TRUE)` -- they propagate Phipson 2016 robust hyperparameter estimation through separate EB steps
- For complex designs, use design matrices without intercept (~ 0 + group)
- `glmTreat()` tests for log fold changes above a threshold (H0: |LFC| <= threshold), preferred over post-hoc LFC filtering
- edgeR v4 introduced bias-corrected APL with `glmQLFit(legacy=FALSE)` as the default -- pre-2025 results are reproducible with `legacy=TRUE`, but the v4 defaults are improved
- `normLibSizes()` is the v4 canonical name (replaces `calcNormFactors()`); `method='TMM'` remains the documented default, with `method='TMMwsp'` as the preferred alternative for sparse / single-cell-pseudobulk data
- For transcript-level DE, use `catchSalmon()` / `catchKallisto()` (Baldoni 2024) instead of standard tximport
- For prokaryotic data, TMM normalization assumes most genes are not DE -- verify this holds

## Related Skills

- deseq2-basics - Alternative DE analysis with DESeq2
- de-visualization - MA plots, volcano plots, heatmaps
- de-results - Extract and export significant genes
- pathway-analysis/go-enrichment - Functional analysis of DE gene lists
