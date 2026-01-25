# DESeq2 Basics - Usage Guide

## Overview

This skill covers differential expression analysis using DESeq2, the most widely-used Bioconductor package for analyzing RNA-seq count data. DESeq2 uses a negative binomial model to test for differential expression between experimental conditions.

## Prerequisites

```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')
BiocManager::install(c('DESeq2', 'apeglm'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Run DESeq2 on my count matrix with treated vs control comparison"
- "Analyze differential expression controlling for batch effects"
- "Get significantly differentially expressed genes with padj < 0.05"

## Example Prompts

### Basic Analysis
> "Create a DESeqDataSet from my count matrix and sample metadata"

> "Run the standard DESeq2 workflow on this RNA-seq data"

> "Apply log fold change shrinkage to my DESeq2 results"

### Design Formulas
> "Set up DESeq2 with batch correction"

> "Create an interaction model for genotype and treatment"

> "Compare multiple conditions against a control"

### Results
> "Extract genes with adjusted p-value < 0.05 and |log2FC| > 1"

> "Order results by significance"

> "Export normalized counts and DE results"

## What the Agent Will Do

1. Create DESeqDataSet from counts and metadata
2. Pre-filter low-count genes
3. Set reference level for comparisons
4. Run DESeq() pipeline
5. Apply lfcShrink() for fold change shrinkage
6. Extract and filter significant results

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| Count matrix | Integer matrix | Genes (rows) x Samples (columns) |
| Sample metadata | Data frame | Rownames matching count column names |
| Design formula | R formula | Variables from metadata (~condition) |

## Tips

- Always use biological replicates (minimum 3 per condition)
- Pre-filter genes with very low counts before analysis
- Use `lfcShrink()` with type='apeglm' for better fold change estimates
- Use `vst()` instead of `rlog()` for large datasets (>100 samples)
- Check `resultsNames(dds)` to see available coefficients for results()
