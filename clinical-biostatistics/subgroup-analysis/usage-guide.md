# Subgroup Analysis - Usage Guide

## Overview

Performs stratified and subgroup analyses for clinical trial data. Covers Mantel-Haenszel pooling across strata, Breslow-Day homogeneity testing, interaction terms in logistic regression, multiple comparisons correction, and forest plot visualization of subgroup effects.

## Prerequisites

```bash
pip install statsmodels scipy numpy pandas matplotlib
```

## Quick Start

Tell your AI agent what you want to do:
- "Test whether the treatment effect differs by age group using an interaction term"
- "Run a Mantel-Haenszel stratified analysis across study sites"
- "Create a forest plot showing treatment effects in each subgroup"
- "Apply multiplicity correction to my subgroup p-values"

## Example Prompts

### Stratified Analysis

> "I have a clinical trial dataset with treatment arm, outcome, and study site. Run a Mantel-Haenszel analysis to get the pooled odds ratio across sites and test whether the OR is consistent across strata with Breslow-Day."

> "Compute the Mantel-Haenszel pooled odds ratio for treatment vs placebo, stratified by disease severity. Include 95% confidence intervals."

### Interaction Testing

> "Fit a logistic regression with treatment, age group, and their interaction to test whether the treatment effect is modified by age. Extract subgroup-specific odds ratios."

> "Test for both multiplicative and additive interaction between treatment and sex. Report RERI for the additive interaction."

### Forest Plots

> "Create a forest plot showing the treatment odds ratio and 95% CI for each pre-specified subgroup: age (<65, 65+), sex, race, baseline severity, and region."

> "Generate a subgroup forest plot with the overall pooled estimate shown as a reference line."

### Regulatory Considerations

> "I have 8 pre-specified subgroups. Apply Holm correction to the subgroup-specific p-values and flag which remain significant after adjustment."

> "Assess whether my subgroup findings meet EMA credibility criteria: pre-specified, significant interaction, consistent across endpoints."

## What the Agent Will Do

1. Load and inspect the clinical dataset for treatment, outcome, and subgroup variables
2. Construct per-stratum 2x2 contingency tables for stratified analysis
3. Compute the Mantel-Haenszel pooled odds ratio with confidence intervals
4. Run the Breslow-Day test for homogeneity of odds ratios across strata
5. Fit a logistic model with interaction terms to formally test effect modification
6. Apply multiplicity correction (Holm or Benjamini-Hochberg) to subgroup p-values
7. Generate a forest plot visualizing subgroup-specific effects with the overall estimate

## Tips

- Always use an interaction term in a single model to test for subgroup effects. Comparing p-values from separate per-subgroup models is statistically invalid because the models have different power.
- Pre-specify subgroups in your statistical analysis plan before unblinding. Post-hoc subgroup analyses are hypothesis-generating only and carry minimal regulatory weight.
- The Breslow-Day test has low power with few strata. A non-significant result does not prove homogeneity. Always supplement with a forest plot for visual assessment.
- For regulatory submissions, apply FWER correction (Holm) to confirmatory subgroup analyses. FDR correction (Benjamini-Hochberg) is appropriate for exploratory screening.
- Report subgroup-specific odds ratios with confidence intervals, not just p-values. The CI width conveys the precision of the estimate in each subgroup.
- When reporting RERI for additive interaction, note that null multiplicative interaction does not imply null additive interaction. The choice depends on whether the research question concerns relative or absolute risk.

## Related Skills

- clinical-biostatistics/categorical-tests - Chi-square and CMH tests used within strata
- clinical-biostatistics/effect-measures - OR computation and forest plots
- clinical-biostatistics/logistic-regression - Interaction terms in regression models
- clinical-biostatistics/trial-reporting - CONSORT-compliant subgroup reporting
- experimental-design/multiple-testing - General multiplicity correction methods
