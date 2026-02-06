# Pleiotropy Detection - Usage Guide

## Overview

Detect and correct for horizontal pleiotropy in Mendelian randomization analyses. Implements MR-PRESSO for outlier removal, MR-Egger regression for directional pleiotropy detection, and Steiger filtering for instrument directionality.

## Prerequisites

```r
install.packages('remotes')
remotes::install_github('MRCIEU/TwoSampleMR')
remotes::install_github('rondolab/MR-PRESSO')

# For MR-RAPS (archived from CRAN March 2025)
install.packages('MendelianRandomization')
```

## Quick Start

Tell your AI agent what you want to do:
- "Run sensitivity analyses on my MR results to check for pleiotropy"
- "Use MR-PRESSO to detect and remove outlier instruments"
- "Check if the MR-Egger intercept is significant for my analysis"

## Example Prompts

### MR-PRESSO

> "Run MR-PRESSO on my harmonized MR data and identify outlier SNPs"

> "Remove pleiotropic outliers and compare the corrected estimate to the original"

### MR-Egger

> "Test for directional pleiotropy using the Egger intercept"

> "What is the I-squared for my MR-Egger analysis -- is the NOME assumption satisfied?"

### Comprehensive Sensitivity

> "Run all standard MR sensitivity analyses: IVW, Egger, weighted median, MR-PRESSO, and leave-one-out"

> "Generate a summary table comparing all MR methods and pleiotropy diagnostics"

### Steiger Filtering

> "Filter my instruments using the Steiger test and re-run MR"

> "Are any of my instruments acting through a reverse causal pathway?"

## What the Agent Will Do

1. Run MR-PRESSO global test for horizontal pleiotropy
2. Identify and remove outlier instruments
3. Test MR-Egger intercept for directional pleiotropy
4. Compute I-squared for NOME assumption assessment
5. Apply Steiger filtering for instrument directionality
6. Compare causal estimates across robust methods
7. Report comprehensive sensitivity analysis summary

## Tips

- **MR-PRESSO sims** - Use at least 5000 distributions for stable p-values; 1000 is a minimum
- **Egger intercept** - Low power with fewer than 10 instruments; non-significant does not mean no pleiotropy
- **I-squared** - Below 0.9 indicates NOME violation; Egger estimate may be biased toward null
- **Method agreement** - Consistent results across IVW, Egger, median, and MR-PRESSO strengthen causal claims
- **MR-RAPS** - The CRAN package was archived March 2025; use the MendelianRandomization package instead
- **STROBE-MR** - Report all sensitivity analyses, not just favorable results
- **Outlier removal** - After removing MR-PRESSO outliers, re-check heterogeneity

## Related Skills

- causal-genomics/mendelian-randomization - Primary MR analysis workflow
- causal-genomics/fine-mapping - Prioritize causal variants at instrument loci
- population-genetics/association-testing - GWAS for instrument selection
