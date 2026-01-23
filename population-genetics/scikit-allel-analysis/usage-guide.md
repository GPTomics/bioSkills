# scikit-allel Analysis - Usage Guide

## Overview

scikit-allel provides Python data structures and algorithms for population genetics analysis. It's ideal for custom analyses, interactive exploration in Jupyter, and integration with other Python tools.

## Prerequisites

```bash
pip install scikit-allel

# Optional for large files
pip install zarr
```

**Note**: scikit-allel is in maintenance mode. For new projects, consider [sgkit](https://github.com/sgkit-dev/sgkit) for long-term support.

## Quick Start

Tell your AI agent what you want to do:
- "Calculate nucleotide diversity from my VCF"
- "Compute allele frequencies per population"
- "Run PCA on my genetic data in Python"
- "Calculate Fst between populations"
- "Analyze haplotype structure"

## Example Prompts

### Loading and Basic Statistics
> "Load my VCF into scikit-allel and calculate allele frequencies"

> "Compute nucleotide diversity (pi) across the genome"

> "Calculate per-site heterozygosity for each sample"

### Population Comparisons
> "Calculate pairwise Fst between my three populations"

> "Run PCA and plot the first two components"

> "Compute Watterson's theta in sliding windows"

### Selection Analysis
> "Calculate Tajima's D in 10kb windows"

> "Compute iHS scores across chromosome 2"

> "Find regions with unusual allele frequency differentiation"

### Large Data Handling
> "Convert my large VCF to Zarr format for efficient access"

> "Calculate statistics on a Zarr-backed dataset"

> "Process my VCF in chunks to avoid memory issues"

## What the Agent Will Do

1. Load VCF data into appropriate array structures
2. Subset data by samples/populations if specified
3. Calculate requested statistics
4. Handle windowing or genome-wide aggregation
5. Generate visualizations if requested
6. Return results as DataFrames or arrays

## Tips

- Use `zarr` backend for VCFs larger than available RAM
- `GenotypeArray` is for diploid data; `HaplotypeArray` for phased haplotypes
- Convert to `AlleleCountsArray` early for faster frequency calculations
- Filter missing data before calculating statistics
- Many functions accept `pos` arrays for windowed calculations

## Data Structures

| Class | Purpose |
|-------|---------|
| `GenotypeArray` | Diploid genotypes (n_var x n_samp x 2) |
| `HaplotypeArray` | Haploid data (n_var x n_hap) |
| `AlleleCountsArray` | Allele counts (n_var x n_alleles) |

## Quick Reference

```python
import allel
import numpy as np

callset = allel.read_vcf('data.vcf.gz')
gt = allel.GenotypeArray(callset['calldata/GT'])
ac = gt.count_alleles()

pi = allel.sequence_diversity(callset['variants/POS'], ac)
print(f'Nucleotide diversity: {pi:.6f}')
```

## Memory Management

For large VCFs, use Zarr:

```python
allel.vcf_to_zarr('large.vcf.gz', 'data.zarr', fields='*')

import zarr
callset = zarr.open('data.zarr', mode='r')
```

## Resources

- [scikit-allel Documentation](https://scikit-allel.readthedocs.io/)
- [scikit-allel GitHub](https://github.com/cggh/scikit-allel)
- [sgkit (successor)](https://github.com/sgkit-dev/sgkit)
