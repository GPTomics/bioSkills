# Selection Statistics - Usage Guide

## Overview

Selection statistics detect signatures of natural selection in genomic data. Different statistics detect different selection types:

| Statistic | Type Detected | Timescale |
|-----------|---------------|-----------|
| Fst | Population differentiation | Any |
| Tajima's D | Neutral departures | Recent |
| iHS | Ongoing sweep | Very recent |
| XP-EHH | Completed sweep | Recent |
| H12/H2H1 | Soft sweeps | Recent |

## When to Use This Skill

- Identifying selection targets
- Comparing selection between populations
- Validating candidate genes
- Understanding evolutionary history

## Installation

```bash
pip install scikit-allel
conda install -c bioconda vcftools
```

## Quick Reference

### Positive Selection (Hard Sweep)

Signs:
- Low Tajima's D (< -2)
- High |iHS| (> 2)
- High Fst between populations
- Reduced diversity (Pi)

### Balancing Selection

Signs:
- High Tajima's D (> 2)
- Elevated heterozygosity
- Old alleles maintained

### Recent Selection

Use haplotype-based methods:
- iHS for ongoing sweeps
- XP-EHH for completed sweeps

### Ancient Selection

Use diversity-based methods:
- Fst for differentiation
- dN/dS for coding regions

## Interpretation Caveats

- Demographic history mimics selection
- Recombination rate affects EHH statistics
- Multiple testing correction needed
- Functional validation recommended

## Resources

- [Selection Tutorial](https://github.com/cggh/scikit-allel/tree/master/docs)
- [vcftools Manual](https://vcftools.github.io/man_latest.html)
