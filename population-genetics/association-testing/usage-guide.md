# Association Testing - Usage Guide

## Overview

GWAS identifies genetic variants associated with traits. PLINK 2.0's `--glm` command provides unified testing for binary (case-control) and quantitative traits using logistic/linear regression with covariate support.

## When to Use This Skill

- Testing SNP associations with disease status
- Quantitative trait mapping
- Candidate gene association studies
- Whole-genome association scans

## Installation

```bash
conda install -c bioconda plink2

# For visualization
pip install pandas matplotlib scipy
```

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
| Genome-wide | 5×10⁻⁸ | Standard GWAS threshold |
| Suggestive | 1×10⁻⁵ | Follow-up candidates |
| Nominal | 0.05 | Not reliable for GWAS |

## Common Issues

### High Genomic Inflation (λ > 1.1)

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
