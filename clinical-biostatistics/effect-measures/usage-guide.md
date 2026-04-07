# Treatment Effect Measures - Usage Guide

## Overview

Computes and interprets treatment effect measures for clinical trial data. Covers odds ratios, risk ratios, and number needed to treat from 2x2 tables and logistic regression models. Includes confidence interval methods, non-collapsibility of ORs, and forest plot visualization for comparing effects across subgroups or studies.

## Prerequisites

```bash
pip install statsmodels numpy pandas matplotlib
```

## Quick Start

Tell your AI agent what you want to do:
- "Compute the odds ratio and risk ratio from my treatment vs placebo 2x2 table with confidence intervals"
- "Extract adjusted odds ratios from a logistic regression model"
- "Calculate the number needed to treat from my clinical trial results"
- "Create a forest plot of odds ratios across subgroups"

## Example Prompts

### Odds Ratios

> "I have a 2x2 table with 45 treated responders, 55 treated non-responders, 30 control responders, and 70 control non-responders. Compute the odds ratio with a 95% confidence interval."

> "Fit a logistic regression predicting response from treatment arm adjusted for age and sex, then extract the adjusted odds ratios."

### Risk Measures

> "Compute both the odds ratio and risk ratio for my clinical trial. The outcome prevalence is 25%, so I want to see how much the OR overestimates the RR."

> "Calculate the NNT for my treatment given an odds ratio of 0.6 and a baseline risk of 15%."

### Visualization

> "Create a forest plot showing the treatment odds ratio separately for males and females, ages under/over 65, and the overall estimate."

> "Plot the subgroup-specific ORs with confidence intervals on a log scale with a reference line at 1.0."

### Interpretation

> "My adjusted and unadjusted ORs differ substantially. Help me determine whether this is due to confounding or non-collapsibility."

> "Compare the NNT at different baseline risk levels for an OR of 0.5 to show stakeholders how the treatment benefit varies by patient risk."

> "What does an odds ratio of 1.52 actually mean in practical terms for my patients? How do I explain it to a clinician?"

## What the Agent Will Do

1. Load the dataset and construct 2x2 contingency tables
2. Compute crude OR and RR with confidence intervals using Table2x2
3. Fit a logistic regression model if covariate adjustment is needed
4. Extract adjusted ORs by exponentiating model coefficients
5. Calculate NNT from absolute risk reduction
6. Generate forest plots for visual comparison across subgroups

## Tips

- Always report confidence intervals alongside point estimates. An OR without a CI is uninterpretable.
- When outcome prevalence exceeds 10%, prefer risk ratios over odds ratios for communicating treatment effects to non-statistical audiences.
- ORs are often reported as percentage changes: `(OR - 1) * 100%`. Be clear this is percentage change in *odds*, not risk -- clinicians often conflate these.
- NNT is only meaningful with a specified baseline risk. The same OR produces very different NNTs depending on how common the outcome is.
- Use a log scale for forest plots so that reciprocal effects (OR 0.5 and OR 2.0) appear equidistant from the null.
- Differences between adjusted and unadjusted ORs may reflect non-collapsibility rather than confounding. This is a mathematical property of the OR, not a bias.
- For small samples, consider profile likelihood confidence intervals instead of Wald intervals, which have poor coverage near boundaries.

## Related Skills

- clinical-biostatistics/logistic-regression - Fit models that produce adjusted ORs
- clinical-biostatistics/categorical-tests - Contingency table tests that complement effect measures
- clinical-biostatistics/subgroup-analysis - Forest plots for subgroup effects
- clinical-biostatistics/trial-reporting - CONSORT-compliant effect reporting
- machine-learning/survival-analysis - Hazard ratios for time-to-event outcomes
