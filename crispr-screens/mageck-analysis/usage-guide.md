# MAGeCK Analysis Usage Guide

## Overview

MAGeCK is the standard tool for analyzing pooled CRISPR screens. It handles count normalization, identifies significantly enriched/depleted genes, and performs pathway analysis.

## Key Concepts

### Negative Selection (Dropout)
Genes where sgRNA depletion indicates essentiality or sensitivity.

### Positive Selection (Enrichment)
Genes where sgRNA enrichment indicates resistance or growth advantage.

### RRA Score
Robust Rank Aggregation score combining sgRNA-level statistics to gene level.

### Beta Score (MLE)
Effect size from maximum likelihood estimation, comparable across genes.

## Workflow Decision

### Use `mageck test` when:
- Simple two-group comparison
- No complex experimental design
- Quick analysis needed

### Use `mageck mle` when:
- Multiple conditions/timepoints
- Need to account for covariates
- Want beta scores for effect sizes

## Normalization Methods

| Method | Use When |
|--------|----------|
| median | Default, most robust |
| total | Even library distribution |
| control | Non-targeting controls available |

## Quality Thresholds

- sgRNA mapping rate: >70%
- Gini index: <0.2 (even distribution)
- Zero-count sgRNAs: <1%
- Replicate correlation: >0.8

## References

- MAGeCK paper: doi:10.1186/s13059-014-0554-4
- MAGeCK-MLE: doi:10.1186/s13059-015-0843-6
- Documentation: https://sourceforge.net/p/mageck/wiki/Home/
