# Linkage Disequilibrium - Usage Guide

## Overview

Linkage disequilibrium (LD) measures non-random association between alleles at different loci. High LD indicates variants are inherited together. LD pruning removes correlated variants for unbiased population structure analysis.

## When to Use This Skill

- Preparing data for PCA or Admixture
- Identifying independent GWAS signals
- Understanding haplotype structure
- Detecting recombination patterns

## Installation

```bash
conda install -c bioconda plink plink2 vcftools
pip install scikit-allel matplotlib
```

## Key Statistics

| Statistic | Range | Interpretation |
|-----------|-------|----------------|
| r² | 0-1 | Correlation squared; 1 = perfect LD |
| D' | 0-1 | Normalized LD; 1 = no recombination |

## Quick Reference

### LD Pruning for PCA

```bash
plink2 --bfile data --indep-pairwise 50 10 0.1 --out prune
plink2 --bfile data --extract prune.prune.in --make-bed --out pruned
```

### Calculate r² Between SNPs

```bash
plink2 --bfile data --r2 --ld-window-kb 500 --out ld
```

### GWAS Clumping

```bash
plink --bfile data --clump results.txt --clump-r2 0.1 --out clumped
```

## Choosing Pruning Thresholds

| Application | r² Threshold | Notes |
|-------------|--------------|-------|
| PCA/Admixture | 0.1 | Strict, independent SNPs |
| GWAS clumping | 0.1-0.2 | Independent signals |
| Polygenic scores | 0.5 | Retain more signal |

## Resources

- [PLINK LD Documentation](https://www.cog-genomics.org/plink/2.0/ld)
- [LD Theory](https://www.nature.com/articles/nrg1123)
