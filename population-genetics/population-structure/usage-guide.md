# Population Structure - Usage Guide

## Overview

Population structure analysis identifies genetic ancestry and stratification using PCA (clustering in continuous space) and ADMIXTURE (discrete ancestry proportions). Essential for GWAS stratification control and ancestry analysis.

## When to Use This Skill

- Assessing population stratification before GWAS
- Identifying sample ancestry
- Detecting outliers or sample swaps
- Understanding genetic diversity
- Generating covariates for association testing

## Installation

```bash
# PLINK 2.0
conda install -c bioconda plink2

# ADMIXTURE
conda install -c bioconda admixture

# Visualization
pip install pandas matplotlib
```

## Quick Start

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
