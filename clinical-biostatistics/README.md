# clinical-biostatistics

## Overview

Statistical analysis methods for clinical trial data, from CDISC data handling through logistic regression, categorical testing, and regulatory-compliant reporting.

**Tool type:** python | **Primary tools:** statsmodels, scipy, tableone, pyreadstat

## Skills

| Skill | Description |
|-------|-------------|
| cdisc-data-handling | Read CDISC SDTM domains, join on USUBJID, aggregate events to subject level |
| logistic-regression | Binary, ordinal, multinomial logistic regression with OR extraction and Firth's method |
| categorical-tests | Chi-square, Fisher's exact, CMH, McNemar's with effect sizes |
| effect-measures | Odds ratios, risk ratios, NNT/NNH, confidence intervals, forest plots |
| subgroup-analysis | Mantel-Haenszel stratification, interaction terms, multiplicity correction |
| trial-reporting | Table 1, ITT/PP populations, multiple imputation, CONSORT/ICH E9 compliance |

## Example Prompts

- "Load my CDISC .xpt files and create a subject-level analysis dataset"
- "Run logistic regression on treatment vs outcome controlling for age and sex"
- "Test association between vaccination status and disease severity with chi-square"
- "Compute adjusted odds ratios with 95% confidence intervals from my clinical trial"
- "Analyze treatment effects across patient subgroups and generate a forest plot"
- "Create a Table 1 of baseline characteristics and handle missing data with multiple imputation"

## Requirements

```bash
pip install statsmodels scipy pingouin tableone pyreadstat pandas numpy matplotlib scikit-learn
```

Optional for rare events / separation:
```bash
pip install firthmodels
```

## Related Skills

- **machine-learning** - Survival analysis and predictive modeling
- **experimental-design** - Power analysis, sample size, multiple testing
- **epidemiological-genomics** - Genomic epidemiology and AMR surveillance
- **workflows** - Clinical trial analysis pipeline
