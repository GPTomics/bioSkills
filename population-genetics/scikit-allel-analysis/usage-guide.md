# scikit-allel Analysis - Usage Guide

## Overview

scikit-allel provides Python data structures and algorithms for population genetics analysis. It's ideal for custom analyses, interactive exploration, and integration with other Python tools.

**Note**: scikit-allel is in maintenance mode. For new projects, consider [sgkit](https://github.com/sgkit-dev/sgkit) for long-term support.

## When to Use This Skill

- Custom population genetics analyses
- Interactive exploration in Jupyter
- Integration with other Python tools
- Need fine-grained control over calculations
- Visualization with matplotlib

## Installation

```bash
pip install scikit-allel

# Optional for large files
pip install zarr
```

## Quick Start

```python
import allel
import numpy as np

callset = allel.read_vcf('data.vcf.gz')
gt = allel.GenotypeArray(callset['calldata/GT'])
ac = gt.count_alleles()

pi = allel.sequence_diversity(callset['variants/POS'], ac)
print(f'Nucleotide diversity: {pi:.6f}')
```

## Data Structures

| Class | Purpose |
|-------|---------|
| `GenotypeArray` | Diploid genotypes (n_var × n_samp × 2) |
| `HaplotypeArray` | Haploid data (n_var × n_hap) |
| `AlleleCountsArray` | Allele counts (n_var × n_alleles) |

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
