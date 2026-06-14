# GWAS Pipeline - Usage Guide

## Overview

This workflow performs genome-wide association studies (GWAS) to identify genetic variants associated with traits or diseases.

## Prerequisites

```bash
conda install -c bioconda plink plink2

# R packages for visualization
install.packages(c('qqman', 'ggplot2'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Run a GWAS on my genotype data"
- "Find variants associated with my phenotype"
- "Perform case-control association testing"

## Example Prompts

### GWAS workflow
> "Run QC and association testing on my VCF"

> "Create Manhattan and QQ plots for my GWAS"

> "Adjust for population structure using PCA"

### Analysis options
> "Run GWAS for a quantitative trait"

> "Include age and sex as covariates"

> "Extract genome-wide significant hits"

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| Genotypes | VCF or PLINK | SNP genotype data |
| Phenotypes | Text file | Case/control or quantitative |
| Covariates | Text file | Age, sex, PCs (optional) |

## What the Workflow Does

1. **QC Filtering** - Remove poor quality samples/variants in order (variant missingness before sample missingness), with controls-only HWE and KING relatedness pruning
2. **Phasing and Imputation** - Align to a reference panel, phase, impute to dosages, and filter by R2 (owned by the phasing-imputation skills; impute cases and controls together)
3. **LD Pruning** - Get independent variants for PCA after excluding long-range-LD regions
4. **PCA** - Calculate population structure covariates
5. **Association** - Test variant-phenotype associations on dosages (PC-covariate GLM for unrelated samples, a linear mixed model for related/structured ones)
6. **Visualization** - Manhattan and QQ plots

## Case-Control vs Quantitative

| Feature | Case-Control | Quantitative |
|---------|--------------|--------------|
| Phenotype | 1=control, 2=case | Continuous value |
| Model | Logistic regression | Linear regression |
| Output | Odds ratio | Beta coefficient |

## Tips

- **Sample size**: Need thousands of samples for common variants
- **Lambda**: An elevated lambda is expected under a polygenic trait and is mostly true signal; do not genomic-control on it. Use the LDSC intercept (or attenuation ratio) to tell polygenicity from confounding
- **Mixed models**: PC covariates absorb continuous ancestry but cannot remove relatedness; use a linear mixed model (BOLT-LMM, SAIGE, regenie) with LOCO for related or finely-structured samples, and SPA or Firth at extreme case:control imbalance
- **Rare variants**: Single-variant tests are underpowered at low minor allele count; aggregate by gene with burden/SKAT/SKAT-O (see population-genetics/rare-variant-association)
- **Multiple testing**: Genome-wide threshold is p < 5e-8 (tighter for WGS and rare variants)
- **Replication**: Always validate findings in independent cohort
- **Imputation**: Impute cases and controls together and carry dosages (not hard calls) into association; see the phasing-imputation skills
