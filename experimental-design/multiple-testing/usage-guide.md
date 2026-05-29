# Multiple Testing Correction Usage Guide

## Overview

This guide covers controlling error rates when thousands of features are tested simultaneously in genomics discovery. The central decision is which error rate to control: the false discovery rate (FDR) for discovery, where tolerating a small controlled fraction of false positives buys large power, versus the family-wise error rate (FWER) for confirmatory work with a few pre-specified hypotheses. The guide covers Benjamini-Hochberg and its dependence assumptions, the Benjamini-Yekutieli procedure for arbitrary dependence, Storey's q-value with proportion-of-true-nulls estimation, local FDR, covariate-weighted FDR via IHW, the independent-filtering rule, the false-coverage-rate trap when reporting confidence intervals only for selected hits, and the GWAS genome-wide threshold.

## Prerequisites

```r
# R/Bioconductor
install.packages('BiocManager')
BiocManager::install(c('qvalue', 'IHW'))

# Python
pip install statsmodels scipy
```

## Quick Start

Tell your AI agent what you want to do:
- "Apply Benjamini-Hochberg correction to my 20,000 DESeq2 p-values at FDR 0.05"
- "Should I use BH, BY, or q-value given my test statistics are correlated?"
- "Compute q-values and the proportion of true nulls (pi0)"
- "Use IHW with mean expression as the covariate to gain power"
- "What genome-wide significance threshold should my GWAS use?"

## Example Prompts

### Differential Expression

> "I have p-values from DESeq2 for 20,000 genes. Apply Benjamini-Hochberg, report the number of discoveries at FDR 0.05, and also give q-values with pi0."

> "My test statistics are strongly correlated across genes. Is BH still valid, or should I use BY?"

### Power Recovery

> "I have mean expression for every gene. Use IHW to weight hypotheses and compare the number of discoveries against plain BH."

> "Is it legitimate to filter out low-count genes before testing to gain power?"

### Method Choice and GWAS

> "What is the difference between an adjusted p-value, a q-value, and a local FDR?"

> "Apply the genome-wide significance threshold to my GWAS summary statistics and explain where it comes from."

## What the Agent Will Do

1. Identify the regime (discovery versus confirmatory) and the appropriate error rate.
2. Check the dependence structure to choose between BH and BY.
3. Apply the correction and report discoveries, with pi0 when q-values are used.
4. Apply IHW or assess an independent filter when a power gain is available and legitimate.
5. Flag the statsmodels default-method trap and the false-coverage-rate issue for selected intervals.

## Tips

- FDR is the discovery default; FWER (Bonferroni/Holm) is for small confirmatory panels.
- BH controls FDR under independence or positive dependence; use BY under arbitrary or negative dependence.
- The q-value estimates the proportion of true nulls and is more powerful than BH when most hypotheses are alternatives.
- IHW recovers power only when the covariate is independent of the p-value under the null.
- In Python, statsmodels multipletests defaults to Holm-Sidak, not BH; always pass the method explicitly.
- Filtering before testing helps only if the filter is independent of the test statistic under the null.

## Related Skills

- power-analysis - The FDR target feeds the power calculation
- sample-size - Replicate number depends on the FDR threshold chosen here
- batch-design - Surrogate variables change the effective number of tests
- differential-expression/de-results - Where the padj column is applied to a DE table
- population-genetics/association-testing - GWAS genome-wide significance machinery
- pathway-analysis/go-enrichment - Correcting enrichment p-values
- clinical-biostatistics/multiplicity-graphical - Confirmatory FWER / closed testing for trials
