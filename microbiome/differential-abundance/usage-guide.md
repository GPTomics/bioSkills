# Differential Abundance - Usage Guide

## Overview

Differential abundance testing identifies taxa that differ significantly between experimental groups while accounting for the compositional nature of microbiome data.

## Prerequisites

```bash
# R packages
BiocManager::install(c('ANCOMBC', 'ALDEx2'))
install.packages('Maaslin2')
```

## Quick Start

Tell your AI agent what you want to do:
- "Find differentially abundant taxa between treatment and control"
- "Run ANCOM-BC2 with covariates for my microbiome study"

## Example Prompts

### Simple Comparisons
> "Run ALDEx2 to compare taxa abundance between two groups"

> "Find differentially abundant genera between healthy and diseased samples"

### Complex Designs
> "Run ANCOM-BC2 with age and sex as covariates"

> "Analyze differential abundance with MaAsLin2 for my longitudinal study"

### Filtering and Interpretation
> "Filter taxa present in fewer than 10% of samples before testing"

> "Which taxa have effect size greater than 1 and q-value below 0.05?"

## What the Agent Will Do

1. Filter low-abundance and rare taxa
2. Select appropriate method based on study design
3. Set up statistical model with covariates if needed
4. Run differential abundance test
5. Apply FDR correction
6. Filter by effect size and significance
7. Generate results table and visualizations

## Tips

- ALDEx2 is best for simple two-group comparisons
- ANCOM-BC2 handles complex designs with covariates
- MaAsLin2 is best for longitudinal and mixed effects
- Always use FDR correction (q-value < 0.05)
- Consider effect size, not just p-value

## Methods to Avoid

- Simple t-test (ignores compositionality)
- DESeq2/edgeR alone (designed for RNA-seq)
- LEfSe (outdated, no FDR control)

## Method Comparison

| Method | Best For | Handles Covariates |
|--------|----------|-------------------|
| ALDEx2 | Two-group comparisons | Limited |
| ANCOM-BC2 | Complex designs | Yes |
| MaAsLin2 | Longitudinal, mixed effects | Yes |

## Filtering Recommendations

Before testing:
- Remove taxa in <10% of samples
- Remove taxa with <0.1% mean abundance
- Remove samples with <1000 reads

## Effect Size Interpretation

| Metric | Meaningful Threshold |
|--------|---------------------|
| ALDEx2 effect | >1 or <-1 |
| Log2 fold change | >1 (2-fold difference) |
