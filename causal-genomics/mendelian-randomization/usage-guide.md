# Mendelian Randomization - Usage Guide

## Overview

Estimate causal effects between exposures and outcomes using genetic variants as instrumental variables. TwoSampleMR implements IVW, MR-Egger, weighted median, and MR-PRESSO methods for robust causal inference from GWAS summary statistics.

## Prerequisites

```r
install.packages('remotes')
remotes::install_github('MRCIEU/TwoSampleMR')

# For OpenGWAS access (see ieugwasr README for current auth instructions)
remotes::install_github('MRCIEU/ieugwasr')
```

## Quick Start

Tell your AI agent what you want to do:
- "Test whether BMI causally affects type 2 diabetes using GWAS summary statistics"
- "Run a two-sample MR analysis with my exposure and outcome GWAS files"
- "Check if cholesterol levels have a causal effect on coronary heart disease"

## Example Prompts

### Basic MR Analysis

> "I have GWAS summary stats for BMI and T2D -- run a Mendelian randomization analysis"

> "Use TwoSampleMR to test if LDL cholesterol causally affects stroke risk"

### Sensitivity Analyses

> "Run MR-Egger, weighted median, and IVW on my harmonized data and compare the results"

> "Check if my MR results are robust using leave-one-out analysis"

### Bidirectional MR

> "Test whether the causal direction is from exposure to outcome using the Steiger test"

> "Run MR in both directions between depression and BMI"

### Local Data

> "I have local GWAS files for exposure and outcome -- set up the MR pipeline without using OpenGWAS"

## What the Agent Will Do

1. Load exposure GWAS data and select genome-wide significant instruments
2. Clump instruments to remove LD (r2 < 0.001, 10 Mb window)
3. Extract matching SNPs from outcome GWAS
4. Harmonize alleles between exposure and outcome
5. Run IVW, MR-Egger, weighted median, and weighted mode methods
6. Assess heterogeneity (Cochran's Q) and pleiotropy (Egger intercept)
7. Generate scatter, forest, funnel, and leave-one-out plots
8. Report causal estimate with confidence intervals

## Tips

- **F-statistic** - Each instrument should have F > 10 to avoid weak instrument bias
- **Method agreement** - Consistent results across IVW, Egger, and median strengthen causal claims
- **Egger intercept** - A significant intercept (p < 0.05) indicates directional pleiotropy
- **OpenGWAS auth** - Auth requirements change; prefer local GWAS files for reproducibility
- **Palindromic SNPs** - Use action = 2 in harmonise_data to handle strand ambiguity
- **Steiger filtering** - Always run directionality_test to confirm instruments act on exposure first
- **Reporting** - Follow STROBE-MR guidelines for transparent reporting of MR studies

## Related Skills

- causal-genomics/pleiotropy-detection - Sensitivity analyses for pleiotropy
- causal-genomics/colocalization-analysis - Confirm shared causal variants at loci
- causal-genomics/fine-mapping - Prioritize causal variants at instrument loci
- population-genetics/association-testing - GWAS for generating summary statistics
