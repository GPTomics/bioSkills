# Association Testing - Usage Guide

## Overview

GWAS identifies genetic variants associated with traits using regression models. PLINK 2.0's `--glm` command provides unified testing for binary (case-control) and quantitative traits with covariate support.

## Prerequisites

```bash
conda install -c bioconda plink2

# For visualization
pip install pandas matplotlib scipy
```

## Quick Start

Tell your AI agent what you want to do:
- "Run a GWAS for my case-control phenotype"
- "Test association with a quantitative trait"
- "Perform association testing with population covariates"
- "Find genome-wide significant hits"
- "Generate a Manhattan plot from my GWAS results"

## Example Prompts

### Basic Association
> "Run a genome-wide association study on my case-control data"

> "Test SNP associations with my quantitative phenotype"

> "Perform logistic regression GWAS for disease status"

### With Covariates
> "Run GWAS including age, sex, and the first 5 PCs as covariates"

> "Test associations while controlling for population stratification"

> "Perform association testing with a custom covariate file"

### Results Analysis
> "Extract all genome-wide significant variants from my GWAS"

> "Calculate genomic inflation factor for my results"

> "Find the top 100 hits and annotate them with gene names"

> "Create Manhattan and QQ plots from my association results"

## What the Agent Will Do

1. Verify input data quality (post-QC PLINK files)
2. Check phenotype file format and distribution
3. Generate PCs for stratification control if needed
4. Run association testing with appropriate model
5. Filter results by significance threshold
6. Calculate genomic inflation (lambda)
7. Generate visualization if requested

## Tips

- Always include population PCs as covariates to control stratification
- Start with `--glm hide-covar` to simplify output
- Genome-wide significance is 5e-8; suggestive is 1e-5
- Lambda > 1.1 suggests residual stratification or relatedness
- For family data, use mixed models (GCTA, BOLT-LMM) instead

## Standard GWAS Workflow

### 1. Quality Control

```bash
plink2 --bfile raw \
    --maf 0.01 --geno 0.05 --mind 0.05 --hwe 1e-6 \
    --make-bed --out qc
```

### 2. Calculate PCs for Stratification

```bash
plink2 --bfile qc --pca 10 --out pca
```

### 3. Run Association

```bash
plink2 --bfile qc \
    --pheno phenotypes.txt \
    --covar pca.eigenvec \
    --covar-name PC1-PC5 \
    --glm hide-covar \
    --out gwas
```

### 4. Identify Significant Hits

```bash
awk '$13 < 5e-8' gwas.PHENO1.glm.* > significant.txt
```

## Significance Thresholds

| Level | P-value | Use |
|-------|---------|-----|
| Genome-wide | 5e-8 | Standard GWAS threshold |
| Suggestive | 1e-5 | Follow-up candidates |
| Nominal | 0.05 | Not reliable for GWAS |

## Common Issues

### High Genomic Inflation (lambda > 1.1)

- Add more PCs as covariates
- Check for cryptic relatedness
- Consider mixed models (GCTA, BOLT-LMM)

### No Significant Results

- Check phenotype file format
- Verify sample sizes
- May be underpowered

### Separation Issues (Logistic)

- Firth regression automatically applied
- Check for very rare variants

## Resources

- [PLINK 2.0 GLM Documentation](https://www.cog-genomics.org/plink/2.0/assoc)
- [GWAS Tutorial](https://github.com/MareesAT/GWA_tutorial)
