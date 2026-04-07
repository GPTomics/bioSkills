# Clinical Trial Analysis Pipeline - Usage Guide

## Overview

End-to-end workflow for analyzing clinical trial data. Takes CDISC SDTM domain files through data preparation, primary analysis (logistic regression), categorical tests, subgroup analysis, missing data handling, and regulatory-compliant reporting.

## Prerequisites

```bash
pip install statsmodels scipy tableone pyreadstat pandas numpy matplotlib scikit-learn
```

Optional for rare events:
```bash
pip install firthmodels
```

## Quick Start

Tell your AI agent what you want to do:
- "Analyze my clinical trial data from start to finish"
- "I have CDISC .xpt files for a vaccine trial -- run the full analysis"
- "Perform a complete statistical analysis of treatment vs placebo on adverse events"
- "Run logistic regression, subgroup analysis, and generate Table 1 for my trial data"

## Example Prompts

### Full Pipeline

> "I have DM, AE, and EX domain files from a BCG vaccination trial. Create a subject-level dataset, run logistic regression on COVID-19 severity controlling for age and patient interaction, test associations with chi-square, analyze subgroups by patient load, and generate a Table 1."

> "Analyze my clinical trial data end to end. The primary endpoint is a binary adverse event. I need adjusted odds ratios, subgroup forest plots, and a sensitivity analysis with multiple imputation for missing baseline data."

### Data Preparation

> "Load my CDISC SDTM files and create a subject-level analysis dataset. Merge demographics with adverse events, aggregate to one row per subject with maximum severity, and code the treatment variable."

### Primary Analysis

> "Run logistic regression on my prepared clinical dataset with treatment as the primary predictor, adjusting for age and sex. Extract odds ratios with 95% confidence intervals."

### Subgroup Analysis

> "Test whether the treatment effect varies across patient subgroups using interaction terms. Generate a forest plot showing subgroup-specific odds ratios with multiplicity-adjusted p-values."

### Reporting

> "Generate a Table 1 of baseline characteristics by treatment arm, run multiple imputation for missing data, and summarize results following CONSORT guidelines."

## What the Agent Will Do

1. Load CDISC domain files (.xpt or .csv) and identify key variables (USUBJID, ARM, outcome)
2. Aggregate event-level data to subject level and merge domains into a single analysis dataset
3. Generate Table 1 baseline characteristics with standardized mean differences
4. Fit logistic regression with explicit reference category and covariate adjustment
5. Run chi-square or Fisher's exact tests on key categorical associations
6. Test subgroup interactions and compute stratum-specific odds ratios
7. Generate a forest plot visualizing subgroup effects
8. Perform multiple imputation sensitivity analysis with Rubin's rules
9. Summarize results with effect estimates, confidence intervals, and CONSORT checklist

## Tips

- Always set an explicit reference category in logistic regression (e.g., Placebo) to avoid incorrect effect direction
- Check expected cell counts before using chi-square; switch to Fisher's exact when any cell < 5
- Use interaction terms to test subgroup effects, not separate per-subgroup models
- Set `sample_posterior=True` in IterativeImputer -- without it, multiple imputation does not properly capture uncertainty
- Report standardized mean differences (SMD > 0.1) rather than p-values for Table 1 balance assessment
- Pre-specify subgroups before unblinding; post-hoc subgroup results are hypothesis-generating only
- When outcome prevalence exceeds 10%, note that odds ratios overestimate risk ratios
- Run the pipeline in order. Data preparation issues (duplicate subjects, unmapped columns) cascade into every downstream step, so validate the analysis dataset before fitting any models

## Related Skills

- clinical-biostatistics/cdisc-data-handling - Detailed CDISC domain handling
- clinical-biostatistics/logistic-regression - Advanced regression methods
- clinical-biostatistics/categorical-tests - Statistical tests for categorical data
- clinical-biostatistics/effect-measures - Effect size computation and interpretation
- clinical-biostatistics/subgroup-analysis - Stratified and interaction analyses
- clinical-biostatistics/trial-reporting - Regulatory standards and missing data
- reporting/rmarkdown-reports - Formatted report generation
