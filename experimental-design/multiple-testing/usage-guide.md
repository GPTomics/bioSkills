# Multiple Testing Correction Usage Guide

## Overview

This guide covers applying multiple testing corrections including FDR, Bonferroni, and q-value methods.

## Prerequisites

```r
# R/Bioconductor
install.packages('BiocManager')
BiocManager::install('qvalue')

# Python
pip install statsmodels scipy
```

## Quick Start

Tell your AI agent what you want to do:
- "Apply FDR correction to my differential expression p-values"
- "Which multiple testing correction should I use for my GWAS results?"
- "Calculate q-values for my DE results"
- "Filter my results at FDR < 0.05"

## Example Prompts

### Differential Expression

> "I have p-values from DESeq2 for 20,000 genes. Apply Benjamini-Hochberg correction and filter at FDR 0.05"

> "Compare the number of significant genes using Bonferroni vs BH correction"

### Method Selection

> "I'm doing an exploratory analysis. Should I use FDR 0.05 or 0.10?"

> "What's the difference between adjusted p-value and q-value?"

### GWAS

> "Apply genome-wide significance threshold to my GWAS results"

## What the Agent Will Do

1. Identify the analysis context
2. Select appropriate correction method
3. Apply correction to p-values
4. Report number of significant results
5. Explain interpretation of corrected values

## Tips

- BH/FDR is standard for most genomics analyses
- Bonferroni is appropriate for small, targeted gene sets
- q-value provides more power than BH when many true positives exist
- FDR 0.05 means 5% of significant calls are expected to be false
- For exploratory work, FDR 0.10 is acceptable
- GWAS uses genome-wide threshold of 5e-8 (Bonferroni for ~1M tests)
