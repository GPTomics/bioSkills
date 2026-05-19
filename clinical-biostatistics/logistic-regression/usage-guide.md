# Logistic Regression - Usage Guide

## Overview

Performs logistic regression for binary and ordinal clinical trial endpoints with explicit declaration of marginal vs conditional estimand per FDA 2023 Final Guidance "Adjusting for Covariates in Randomized Clinical Trials." Covers covariate adjustment, g-computation/standardisation for marginal effects, modified Poisson regression for direct RR estimation when prevalence > 10%, Brant test for proportional odds, Firth penalty for separation/rare events, and Hauck-Donner detection.

## Prerequisites

```bash
pip install statsmodels scipy numpy pandas scikit-learn firthmodels
```

R is recommended for production marginal effects:

```r
install.packages(c('marginaleffects', 'RobinCar', 'riskCommunicator', 'brant', 'VGAM', 'MASS'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Logistic regression with treatment + covariates; report BOTH conditional OR and marginal RD via g-computation per FDA 2023"
- "Modified Poisson with HC3 sandwich SE to estimate RR directly when prevalence > 10%"
- "Brant test for proportional-odds assumption in ordinal model; partial PO model if violated"
- "Firth penalised logistic regression with penalised LR test (NOT Wald) when separation detected"
- "Detect Hauck-Donner pathology and switch to profile-likelihood CI"

## Example Prompts

### Standard adjusted analysis

> "Fit a logistic regression for response from treatment, age, sex, baseline_severity. Use Placebo as the reference category explicitly (do NOT rely on alphabetical default which gives 'Active' as reference)."

> "Extract conditional ORs for each predictor and the marginal RD via g-computation. Report HC3 sandwich SE for the marginal estimand per FDA 2023."

### Marginal vs conditional (FDA 2023)

> "Per FDA 2023 covariate adjustment guidance, my primary estimand is marginal RD. The model-fit OR is the CONDITIONAL log-OR (a different parameter). Compute marginal RD via g-computation in marginaleffects (R) or manual standardisation (Python)."

> "Explain to my reviewer why the conditional OR and marginal RD differ -- cite Permutt 2020 non-collapsibility."

### Common outcomes -- direct RR

> "Outcome prevalence is 32%, so OR will substantially overstate RR. Switch to modified Poisson regression with cov_type='HC1' (Zou 2004) to estimate RR directly."

### Ordinal outcomes

> "Fit a proportional-odds cumulative-logit model for severity (mild/moderate/severe). Then run the Brant test (R brant::brant) to verify PO holds per predictor. If PO fails for one predictor, fit a partial PO model in VGAM."

### Separation and Firth

> "My binary outcome has only 3% events -- separation likely. Detect with firthmodels.detect_separation. Apply Firth penalty; use penalised LR test (firth.pvalues_lrt_) NOT Wald p-values for inference."

### Hauck-Donner

> "Small cell count; coefficient is large but Wald p > 0.05. Hauck-Donner suspected. Switch to profile-likelihood CI (R MASS::confint.glm); detect with VGAM::hdeff."

### Stratified randomisation

> "My trial used stratified randomisation by site and baseline severity. Include both as model covariates per Kahan-Morris 2012 -- ignoring inflates Type-I by up to 30%."

## What the Agent Will Do

1. Verify the binary/ordinal outcome and explicitly set reference category for treatment
2. Include pre-specified covariates from the SAP (including all randomisation strata)
3. Fit logistic regression and extract conditional ORs with CIs
4. Compute marginal RD via g-computation per FDA 2023 (HC3 sandwich SE)
5. For prevalence > 10%, fit modified Poisson with sandwich SE for direct RR
6. For ordinal outcomes, run Brant test for PO; fall back to partial PO or multinomial as needed
7. Detect separation/rare events and apply Firth penalty with penalised LR inference
8. Run model diagnostics: ROC-AUC, calibration plot (primary), Hosmer-Lemeshow (supplementary), pseudo R-squared

## Tips

- **Reference category is the single most common silent bug.** statsmodels default is alphabetical -- 'Active' sorts before 'Placebo', silently flipping OR direction. Always pass `C(ARM, Treatment(reference="Placebo"))`.
- **Per FDA 2023, the primary estimand is marginal** (RD via g-computation, HC3 SE). The conditional OR from logistic is a DIFFERENT parameter due to non-collapsibility (Permutt 2020).
- **Tsiatis 2008 robustness:** under randomisation Z ⊥ X, g-computation marginal estimator is consistent for marginal ATE even if outcome model is misspecified. Efficiency depends on model quality; consistency does not.
- **Modified Poisson (Zou 2004) requires sandwich SE** (`cov_type='HC1'` or `'HC3'`). Without it, SEs are wrong because binary data are NOT Poisson.
- **Brant test localises PO violations per coefficient.** If only 1-2 predictors violate PO, fit a partial PO model in VGAM rather than abandoning the model.
- **`OrderedModel` requires no intercept** -- threshold parameters replace it.
- **Firth's `pvalues_` is Wald and is liberal under penalty.** Use `pvalues_lrt_` (firthmodels >=0.3) or compute penalised LR test manually (Heinze-Schemper 2002).
- **Hauck-Donner detection:** large coefficient with small Wald chi-square indicates the test statistic is non-monotone near the boundary. Use profile likelihood or LR test (Yee 2022).
- **Stepwise covariate selection inflates Type-I.** Pre-specify covariates in SAP based on clinical knowledge.
- **Don't adjust for mediators** (post-treatment variables on causal path); attenuates effect toward null. Use causal DAG to distinguish confounder vs mediator.
- **`sm.Logit` provides `pred_table()` and `prsquared`** that `sm.GLM(family=Binomial)` does not. Use Logit unless changing link functions.
- **EPV >= 10** (events per variable; Peduzzi 1996). Below this, switch to Firth or reduce covariates.
- **H-L is supplementary, not primary.** Use calibration plot as primary diagnostic. H-L has low power n<200 and oversensitivity n>2000.

## Related Skills

- clinical-biostatistics/cdisc-data-handling - Prepare ADaM/SDTM analysis datasets
- clinical-biostatistics/effect-measures - Marginal-vs-conditional in depth; CI methods
- clinical-biostatistics/categorical-tests - Unadjusted alternatives
- clinical-biostatistics/subgroup-analysis - Interaction terms and HTE
- clinical-biostatistics/survival-analysis - Cox regression analogue for TTE
- clinical-biostatistics/multiplicity-graphical - Multiple endpoints from one logistic model
- clinical-biostatistics/trial-reporting - CONSORT 2025 + ICH E9(R1) reporting
