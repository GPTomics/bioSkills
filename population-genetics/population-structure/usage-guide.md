# Population Structure - Usage Guide

## Overview

Population structure analysis identifies genetic ancestry and stratification using PCA (continuous clustering) and ADMIXTURE (discrete ancestry proportions). Essential for GWAS stratification control and ancestry inference.

## Prerequisites

```bash
# PLINK 2.0
conda install -c bioconda plink2

# ADMIXTURE
conda install -c bioconda admixture

# FlashPCA2 (for large datasets)
conda install -c bioconda flashpca

# Visualization
pip install pandas matplotlib
```

## Quick Start

Tell your AI agent what you want to do:
- "Run PCA on my genetic data"
- "Run fast PCA on biobank-scale data with FlashPCA2"
- "Estimate ancestry proportions with ADMIXTURE"
- "Check for population stratification before GWAS"
- "Plot PC1 vs PC2 colored by population"

## Example Prompts

### PCA Analysis
> "Calculate the first 10 principal components from my PLINK data"

> "Run PCA and identify any outlier samples"

> "Generate a PCA plot colored by self-reported ancestry"

> "Run FlashPCA2 on my large dataset (100k+ samples)"

### Admixture Analysis
> "Run ADMIXTURE for K=2 through K=6 and find the best K"

> "Estimate ancestry proportions assuming 3 ancestral populations"

> "Create a stacked bar plot of admixture proportions"

### Combined Analysis
> "Perform full population structure analysis with LD pruning, PCA, and ADMIXTURE"

> "Check my data for population stratification and generate covariates for GWAS"

> "Compare PCA clustering with ADMIXTURE assignments"

## What the Agent Will Do

1. LD prune the data for unbiased estimation
2. Run PCA to calculate principal components
3. Run ADMIXTURE across multiple K values if requested
4. Identify optimal K using cross-validation error
5. Generate visualizations (PCA scatter, ADMIXTURE bar plots)
6. Flag outlier samples if detected

## Tips

- Always LD prune before PCA/ADMIXTURE (r2 < 0.1)
- PC1-2 usually capture the largest population splits
- Choose K with lowest cross-validation error
- Outliers may indicate sample swaps, contamination, or unique ancestry
- Remove related individuals (IBD > 0.125) before analysis
- Use FlashPCA2 for biobank-scale data (100k+ samples) for better performance

## Quick Reference

### PCA Only

```bash
# Simple PCA
plink2 --bfile data --pca 10 --out pca
```

### Full Pipeline

```bash
# 1. LD prune
plink2 --bfile data --indep-pairwise 50 10 0.1 --out prune
plink2 --bfile data --extract prune.prune.in --make-bed --out pruned

# 2. PCA
plink2 --bfile pruned --pca 10 --out pca

# 3. Admixture
for K in 2 3 4 5; do
    admixture --cv pruned.bed $K
done
```

## Interpreting Results

### PCA

- **PC1**: Usually largest population split
- **Clusters**: Groups with similar ancestry
- **Outliers**: Sample swaps, contamination, or unique ancestry

### Admixture

- **K**: Number of ancestral populations
- **Q values**: Proportion of each ancestry
- **CV error**: Lower is better for choosing K

## Common Issues

### PCA Shows No Structure

- May be homogeneous population
- Try more PCs
- Check for batch effects

### Admixture Won't Converge

- LD prune more aggressively
- Remove closely related individuals
- Increase iterations

## Resources

- [ADMIXTURE Manual](https://dalexander.github.io/admixture/admixture-manual.pdf)
- [PLINK PCA Documentation](https://www.cog-genomics.org/plink/2.0/strat)
