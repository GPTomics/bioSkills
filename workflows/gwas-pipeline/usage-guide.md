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

1. **QC Filtering** - Remove poor quality samples/variants
2. **LD Pruning** - Get independent variants for PCA
3. **PCA** - Calculate population structure covariates
4. **Association** - Test variant-phenotype associations
5. **Visualization** - Manhattan and QQ plots

## Case-Control vs Quantitative

| Feature | Case-Control | Quantitative |
|---------|--------------|--------------|
| Phenotype | 1=control, 2=case | Continuous value |
| Model | Logistic regression | Linear regression |
| Output | Odds ratio | Beta coefficient |

## Tips

- **Sample size**: Need thousands of samples for common variants
- **Lambda**: Should be ~1.0; high values indicate stratification
- **Multiple testing**: Genome-wide threshold is p < 5e-8
- **Replication**: Always validate findings in independent cohort
