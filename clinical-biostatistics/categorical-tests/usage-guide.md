# Categorical Association Tests - Usage Guide

## Overview

Tests associations between categorical variables in clinical trial data. Covers Pearson chi-square, Fisher's exact, Cochran-Mantel-Haenszel stratified tests, McNemar's test for paired data, effect size computation (phi, Cramer's V), and post-hoc pairwise comparisons with multiple testing correction.

## Prerequisites

```bash
pip install scipy statsmodels pingouin pandas numpy
```

## Quick Start

Tell your AI agent what you want to do:
- "Test whether treatment response differs between Active and Placebo arms using a chi-square test"
- "Run Fisher's exact test on my 2x2 table because expected counts are too small"
- "Perform a Cochran-Mantel-Haenszel test controlling for study site"
- "Compute Cramer's V to measure effect size for my treatment-outcome table"

## Example Prompts

### Chi-Square Analysis

> "I have a clinical dataset with treatment arm and adverse event category. Test whether the distribution of adverse events differs between treatment groups."

> "Run a chi-square test of independence between treatment and response (responder/non-responder). Check expected cell counts and switch to Fisher's exact if any are below 5."

### Fisher's Exact

> "My 2x2 table has small expected counts. Run Fisher's exact test and report the odds ratio with a two-sided p-value."

> "Compare the rate of serious adverse events between arms using Fisher's exact test. I only have 30 subjects per group."

### Stratified Analysis

> "Test whether treatment is associated with response while controlling for study site using the Cochran-Mantel-Haenszel test. Also check if the odds ratio is consistent across sites."

> "I need to adjust my chi-square analysis for a confounding variable. Show me how to use stratified tables."

### Effect Sizes

> "Compute Cramer's V for my treatment-outcome contingency table and interpret the effect size."

> "Run post-hoc pairwise chi-square tests across all outcome categories and adjust p-values using the Holm method."

## What the Agent Will Do

1. Load the dataset and verify categorical variable types
2. Construct a contingency table and inspect expected cell counts
3. Select the appropriate test (chi-square, Fisher's exact, or CMH) based on sample size and expected counts
4. Compute the test statistic and p-value
5. Calculate the relevant effect size (phi or Cramer's V)
6. Run post-hoc pairwise comparisons with multiple testing correction if more than two outcome categories are present

## Tips

- Always check expected cell counts before interpreting a chi-square test. If any expected count is below 5, use Fisher's exact test instead.
- Use `correction=False` in `chi2_contingency()` for standard Pearson chi-square. The default Yates' correction is overly conservative.
- For stratified analyses, verify homogeneity of odds ratios across strata with the Breslow-Day test before interpreting the pooled MH estimate.
- Report effect sizes (Cramer's V) alongside p-values. A statistically significant result with a tiny effect size may not be clinically meaningful.
- Use FWER corrections (Bonferroni, Holm) for confirmatory analyses and FDR corrections for exploratory analyses.
- McNemar's test requires paired data -- it is not appropriate for independent groups.

## Related Skills

- clinical-biostatistics/effect-measures - Detailed OR and risk ratio analysis
- clinical-biostatistics/logistic-regression - Regression-based alternative to contingency tests
- clinical-biostatistics/subgroup-analysis - Stratified analysis with interaction terms
- experimental-design/multiple-testing - General multiple testing correction methods
