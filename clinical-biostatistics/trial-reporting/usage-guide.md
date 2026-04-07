# Trial Reporting - Usage Guide

## Overview

Prepares statistical tables and reports for clinical trials following regulatory standards. Generates Table 1 baseline characteristics with the tableone package, defines analysis populations (ITT, per-protocol, safety), performs multiple imputation for missing data using Rubin's rules, and structures results per CONSORT and ICH E9 guidelines.

## Prerequisites

```bash
pip install tableone statsmodels scikit-learn pandas numpy
```

## Quick Start

Tell your AI agent what you want to do:
- "Generate a Table 1 comparing baseline characteristics between treatment arms"
- "Define ITT, per-protocol, and safety populations from my clinical dataset"
- "Handle missing outcome data with multiple imputation and pool results with Rubin's rules"
- "Check my analysis report against CONSORT statistical reporting requirements"

## Example Prompts

### Table 1

> "Create a Table 1 from my trial data with age, sex, race, BMI, and disease stage. Group by treatment arm, include p-values and standardized mean differences, and export to Excel."

> "Generate baseline characteristics for my RCT. Flag any variables with SMD greater than 0.1 as potentially imbalanced. Include an overall column."

### Missing Data

> "Assess the pattern of missing data in my outcome variable. Run multiple imputation with 20 datasets, fit a logistic regression on each, and pool the results using Rubin's rules."

> "My trial has 15% missing primary endpoint data. Perform multiple imputation assuming MAR, then run a complete-case sensitivity analysis for comparison."

### Analysis Populations

> "Define ITT, per-protocol, and safety populations from my SDTM data. The ITT set is all randomized, per-protocol excludes major protocol violators, and safety includes anyone who received at least one dose."

> "Create a CONSORT flow diagram showing how many subjects were screened, randomized, allocated to each arm, completed, and included in each analysis population."

### Regulatory Reporting

> "Review my clinical trial analysis against the ICH E9(R1) estimand framework. I need to define the estimand for my primary endpoint accounting for treatment discontinuation."

> "Prepare a statistical methods section following CONSORT guidelines. Include sample size justification, primary analysis method, handling of missing data, and multiplicity adjustment strategy."

## What the Agent Will Do

1. Load the clinical dataset and identify treatment arm, outcome, and baseline variables
2. Generate a Table 1 with descriptive statistics, p-values, and SMDs across arms
3. Define analysis populations (ITT, per-protocol, safety) based on disposition data
4. Assess missing data patterns and proportions across key variables
5. Perform multiple imputation with IterativeImputer using sample_posterior=True
6. Fit the analysis model on each imputed dataset and pool results via Rubin's rules
7. Structure results with effect estimates, confidence intervals, and p-values
8. Verify reporting completeness against CONSORT and ICH E9 checklists

## Tips

- Use SMD rather than p-values to assess baseline balance. An SMD > 0.1 suggests meaningful imbalance regardless of statistical significance. CONSORT 2010 discouraged baseline significance tests.
- Always set `sample_posterior=True` in IterativeImputer for multiple imputation. Without it, imputed datasets are nearly identical and confidence intervals will be artificially narrow.
- Use at least 20 imputations. Older guidance suggested 5-10, but modern recommendations favor m >= 20 for stable estimates, especially when the fraction of missing information is substantial.
- ITT should be the primary analysis population per ICH E9. Per-protocol analyses inflate treatment effects by excluding subjects who may have dropped out due to treatment failure.
- Define the estimand (what is being estimated) before choosing the statistical method. The estimand drives the analysis, not vice versa.
- When missing data exceeds 40% on a key variable, treat results as hypothesis-generating. Multiple imputation under MAR becomes increasingly unreliable with high missingness.
- For co-primary endpoints, pre-specify the multiplicity strategy (Bonferroni, hierarchical, Hochberg) in the statistical analysis plan before unblinding.

## Related Skills

- clinical-biostatistics/cdisc-data-handling - Load and prepare the analysis dataset
- clinical-biostatistics/logistic-regression - Primary analysis models
- clinical-biostatistics/effect-measures - Compute and interpret treatment effects
- clinical-biostatistics/subgroup-analysis - Pre-specified subgroup analyses
- reporting/rmarkdown-reports - Generate formatted statistical reports
- experimental-design/multiple-testing - General multiplicity correction
