# Logistic Regression - Usage Guide

## Overview

Fits logistic regression models for binary and ordinal clinical trial endpoints. Covers treatment effect estimation with odds ratios, covariate adjustment, model diagnostics, and Firth penalized regression for rare events or separation.

## Prerequisites

```bash
pip install statsmodels scipy numpy pandas scikit-learn firthmodels
```

## Quick Start

Tell your AI agent what you want to do:
- "Fit a logistic regression to predict treatment response with age and sex as covariates"
- "Calculate odds ratios with confidence intervals for my clinical trial data"
- "Run Firth's logistic regression because I have separation in my binary outcome"
- "Fit an ordinal logistic model for mild/moderate/severe outcomes"

## Example Prompts

### Binary Outcomes

> "Fit a logistic regression predicting adverse event occurrence from treatment arm, adjusting for age and sex. Use Placebo as the reference group."

> "Model the binary endpoint 'responder' using treatment, baseline severity, and site as predictors. Show me the odds ratios table."

> "Compare the odds of serious adverse events between Active and Placebo arms after adjusting for age and baseline BMI"

### Ordinal Outcomes

> "My outcome is a 4-level severity score. Should I dichotomize it to severe vs not-severe, or use ordinal logistic regression?"

> "Fit an ordinal logistic regression for disease severity (mild/moderate/severe) as a function of treatment and age"

> "Model my 5-level pain score as an ordinal outcome with proportional odds regression"

### Covariate Adjustment

> "Add baseline disease severity, age, and sex as covariates to my treatment effect model. Show me which covariates are significant."

> "I need to adjust for potential confounders in my logistic model. Help me decide which variables to include based on my study design."

### Rare Events

> "My outcome has only 3% prevalence and I'm getting convergence warnings. Apply Firth's penalized logistic regression."

> "Check whether my data has separation issues and if so, fit a Firth model instead of standard logistic regression"

## What the Agent Will Do

1. Load the analysis dataset and verify the binary/ordinal outcome variable
2. Set the reference category for the treatment variable explicitly
3. Fit a logistic regression model with pre-specified covariates
4. Extract odds ratios with 95% confidence intervals and p-values
5. Run model diagnostics (ROC-AUC, Hosmer-Lemeshow, pseudo R-squared)
6. Check for separation if convergence issues arise and switch to Firth's method if needed
7. Report results with clinical interpretation of effect direction and magnitude

## Tips

- Always set the reference level for your treatment variable explicitly. Do not rely on alphabetical ordering -- 'Active' sorts before 'Placebo', which silently reverses the odds ratio direction.
- Use the formula API (`smf.logit`) rather than the matrix API for clinical analyses. It handles categorical variables automatically and produces more readable output.
- Choose covariates from your statistical analysis plan or based on clinical knowledge, not from data-driven stepwise selection.
- If your outcome has fewer than 10 events per predictor variable, consider reducing the number of covariates or using Firth's method to avoid overfitting.
- McFadden's pseudo R-squared values above 0.2 indicate excellent fit -- do not compare directly to OLS R-squared values.
- For ordinal models, do not add an intercept term. The threshold parameters in OrderedModel replace the intercept.
- When reporting odds ratios, always state the reference category and the direction of the effect to avoid misinterpretation.

## Related Skills

- clinical-biostatistics/cdisc-data-handling - Prepare analysis datasets from CDISC domains
- machine-learning/survival-analysis - Time-to-event modeling with Cox regression
