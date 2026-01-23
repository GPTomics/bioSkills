# Selection Statistics - Usage Guide

## Overview

Selection statistics detect signatures of natural selection in genomic data. Different methods detect different selection types and timescales, from recent sweeps (iHS) to ancient differentiation (Fst).

## Prerequisites

```bash
pip install scikit-allel
conda install -c bioconda vcftools
```

## Quick Start

Tell your AI agent what you want to do:
- "Calculate Fst between my two populations"
- "Scan for selection signatures using Tajima's D"
- "Compute iHS to detect ongoing selective sweeps"
- "Find regions under balancing selection"
- "Compare selection pressures between populations"

## Example Prompts

### Diversity Statistics
> "Calculate Tajima's D in 50kb windows across the genome"

> "Compute nucleotide diversity (pi) for each population"

> "Find regions with unusually low diversity suggesting sweeps"

### Population Differentiation
> "Calculate Fst between European and African samples"

> "Find highly differentiated SNPs between cases and controls"

> "Generate a Manhattan plot of Fst values"

### Haplotype-Based Tests
> "Compute iHS scores to detect ongoing selection"

> "Run XP-EHH between my populations to find completed sweeps"

> "Identify haplotypes under positive selection"

### Multi-Statistic Analysis
> "Scan for selection using Tajima's D, Fst, and iHS together"

> "Find regions significant in multiple selection tests"

> "Compare selection signatures across chromosomes"

## What the Agent Will Do

1. Assess data format and phase status
2. Calculate requested statistics genome-wide or in windows
3. Standardize/normalize scores where appropriate
4. Identify outlier regions exceeding thresholds
5. Generate visualizations (Manhattan plots, histograms)
6. Report candidate regions with coordinates

## Tips

- Haplotype-based tests (iHS, XP-EHH) require phased data
- Demographic history can mimic selection signals
- Use multiple statistics to reduce false positives
- Always adjust for recombination rate variation
- Empirical outlier cutoffs (top 1%) are often more reliable than p-values

## Selection Signatures Reference

| Statistic | Type Detected | Timescale |
|-----------|---------------|-----------|
| Fst | Population differentiation | Any |
| Tajima's D | Neutral departures | Recent |
| iHS | Ongoing sweep | Very recent |
| XP-EHH | Completed sweep | Recent |
| H12/H2H1 | Soft sweeps | Recent |

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
