# Linkage Disequilibrium - Usage Guide

## Overview

Linkage disequilibrium (LD) measures non-random association between alleles at different loci. LD pruning removes correlated variants for unbiased population structure analysis; LD clumping identifies independent GWAS signals.

## Prerequisites

```bash
conda install -c bioconda plink plink2 vcftools
pip install scikit-allel matplotlib
```

## Quick Start

Tell your AI agent what you want to do:
- "LD prune my data for PCA analysis"
- "Calculate r-squared between SNPs"
- "Clump my GWAS results to find independent signals"
- "Generate an LD heatmap for a candidate region"
- "Find tag SNPs for my variants of interest"

## Example Prompts

### LD Pruning
> "LD prune my PLINK data with r2 threshold 0.1 for population structure analysis"

> "Create an independent SNP set for ADMIXTURE"

> "Prune my data keeping one variant per 50kb window"

### LD Calculation
> "Calculate r-squared between all pairs of SNPs within 500kb"

> "Compute LD for variants in the HLA region"

> "Generate an LD matrix for my candidate locus"

### GWAS Clumping
> "Clump my GWAS results to identify independent signals"

> "Find lead SNPs for each associated locus"

> "Extract independent hits with r2 < 0.1"

### Visualization
> "Create an LD heatmap for chromosome 6p21"

> "Plot LD decay with distance"

> "Visualize haplotype blocks in my region of interest"

## What the Agent Will Do

1. Determine appropriate LD operation (prune, calculate, or clump)
2. Set window size and r2 threshold based on application
3. Run PLINK or scikit-allel commands
4. Report number of variants before/after pruning
5. Generate visualization if requested

## Tips

- Use r2 < 0.1 for PCA/ADMIXTURE (strict independence)
- Use r2 < 0.2 for GWAS clumping (independent signals)
- Use r2 < 0.5 for polygenic scores (retain more signal)
- LD patterns vary by population; use matched reference panels
- Include centromeric/telomeric regions in LD analysis can cause artifacts

## Key Statistics

| Statistic | Range | Interpretation |
|-----------|-------|----------------|
| r2 | 0-1 | Correlation squared; 1 = perfect LD |
| D' | 0-1 | Normalized LD; 1 = no recombination |

## Quick Reference

### LD Pruning for PCA

```bash
plink2 --bfile data --indep-pairwise 50 10 0.1 --out prune
plink2 --bfile data --extract prune.prune.in --make-bed --out pruned
```

### Calculate r2 Between SNPs

```bash
plink2 --bfile data --r2 --ld-window-kb 500 --out ld
```

### GWAS Clumping

```bash
plink --bfile data --clump results.txt --clump-r2 0.1 --out clumped
```

## Choosing Pruning Thresholds

| Application | r2 Threshold | Notes |
|-------------|--------------|-------|
| PCA/Admixture | 0.1 | Strict, independent SNPs |
| GWAS clumping | 0.1-0.2 | Independent signals |
| Polygenic scores | 0.5 | Retain more signal |

## Resources

- [PLINK LD Documentation](https://www.cog-genomics.org/plink/2.0/ld)
- [LD Theory](https://www.nature.com/articles/nrg1123)
