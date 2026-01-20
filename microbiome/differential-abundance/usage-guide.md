# Differential Abundance Testing

## Overview

Differential abundance testing identifies taxa that differ significantly between experimental groups while accounting for the compositional nature of microbiome data.

## The Compositionality Problem

Microbiome data from sequencing is compositional:
- Abundances are relative (proportions)
- Changes in one taxon affect apparent abundance of others
- Standard statistical tests give false positives

## Recommended Methods

### ALDEx2
- CLR (centered log-ratio) transformation
- Monte Carlo sampling from Dirichlet distribution
- Handles compositionality properly
- Best for: Simple two-group comparisons

### ANCOM-BC
- Bias correction for compositionality
- Handles structural zeros
- Supports covariates
- Best for: Complex designs with covariates

### MaAsLin2
- Multiple normalization options
- Linear models with covariates
- TSS, CLR, or AST normalization
- Best for: Mixed effects, longitudinal data

## Methods to Avoid

- **Simple t-test**: Ignores compositionality
- **DESeq2/edgeR alone**: Designed for RNA-seq, not compositional
- **LEfSe**: Outdated, no FDR control

## Workflow

1. **Filter low-abundance taxa**: Remove rare taxa
2. **Choose appropriate method**: Based on design
3. **Run differential test**: With proper normalization
4. **Apply FDR correction**: BH or similar
5. **Effect size filtering**: Not just p-value

## Filtering Recommendations

Before testing:
- Remove taxa in <10% of samples
- Remove taxa with <0.1% mean abundance
- Remove samples with <1000 reads

## Interpreting Results

### Effect Size
- **ALDEx2 effect**: >1 or <-1 is meaningful
- **Log2 fold change**: >1 is 2-fold difference

### Multiple Testing
- Always use FDR correction (BH)
- q-value < 0.05 is standard threshold

### Biological Significance
- Consider absolute abundance
- Validate with qPCR if possible
- Check consistency across methods
