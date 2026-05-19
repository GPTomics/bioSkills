# Treatment Effect Measures - Usage Guide

## Overview

Computes and interprets treatment effect measures (OR, RR, RD, HR, NNT) for clinical trial data with calibrated confidence intervals (Wilson, Newcombe, Miettinen-Nurminen, MOVER, profile likelihood, Bender NNT) and explicit declaration of whether the estimand is marginal or conditional under ICH E9(R1) and FDA 2023 covariate adjustment guidance.

## Prerequisites

```bash
pip install statsmodels numpy pandas matplotlib
```

R is recommended for production:

```r
install.packages(c('ratesci', 'exact2x2', 'marginaleffects', 'RobinCar', 'riskCommunicator'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Compute OR and RR from my 2x2 with Wald, profile-likelihood, and Miettinen-Nurminen score CIs"
- "Extract marginal RD via g-computation from a covariate-adjusted logistic regression per FDA 2023"
- "Calculate NNT with Bender 2002 convention (NNTB-infinity-NNTH when RD CI crosses zero)"
- "Forest plot of subgroup ORs on log scale with reference line"
- "Modified Poisson with HC1/HC3 sandwich SE to estimate RR directly when prevalence > 10%"

## Example Prompts

### Crude effect measures

> "I have a 2x2: 45 treated responders, 55 treated non-responders, 30 control responders, 70 control non-responders. Compute OR, RR, RD with calibrated CIs (Wald + Newcombe-Wilson hybrid)."

> "Use Miettinen-Nurminen score CI for the RD as the regulatory standard for NI margin assessment."

### Marginal vs conditional (FDA 2023)

> "Per FDA 2023 covariate adjustment guidance, my primary estimand should be marginal RD. Fit logistic adjusted for age + sex + baseline severity, then compute marginal RD via g-computation with HC3 sandwich SE. Report conditional OR as supportive."

> "My adjusted and unadjusted ORs differ -- is it confounding or non-collapsibility? Cite Permutt 2020."

### NNT with proper convention

> "Compute NNT from my trial with treated 30/200 events, control 50/200. Use Bender 2002 method; if RD CI crosses zero, report NNTB-infinity-NNTH convention."

> "Compare NNT at baseline risks 5%, 20%, 50% for an OR of 0.5 to show stakeholders how benefit varies."

### Modified Poisson for RR

> "Outcome prevalence is 25% -- OR will overstate RR. Use modified Poisson regression with HC1 sandwich SE to estimate RR directly (Zou 2004)."

### Forest plots and visualisation

> "Create a forest plot of subgroup ORs (sex, age <65/>=65, race, severity, region) on log scale with reference line at 1.0."

### Hauck-Donner pathology

> "Small cell counts; Wald p is non-significant but coefficient is large. Switch to profile-likelihood CI (R MASS::confint.glm); cite Yee 2022."

## What the Agent Will Do

1. Construct the 2x2 table with proper column ordering (event first per Table2x2)
2. Compute OR, RR, RD using calibrated CI methods (Wilson, Newcombe-Wilson, MN, profile likelihood)
3. For binary outcomes, fit logistic regression and compute BOTH conditional OR AND marginal RD via g-computation per FDA 2023
4. Apply modified Poisson with sandwich SE when prevalence > 10% and RR is the policy quantity
5. Convert effect to NNT with Bender 2002 disjoint-interval convention when RD CI crosses zero
6. Generate forest plots on log scale for OR/RR/HR; linear for RD/RMST

## Tips

- **OR is non-collapsible.** Adjusted vs unadjusted ORs differ even without confounding (Permutt 2020). Don't interpret the difference as confounding.
- **Per FDA 2023, marginal RD via g-computation is the primary estimand** for binary endpoints with covariate adjustment. Conditional OR is a different parameter, not a different estimate.
- **Wald CI has well-documented coverage failures** (Brown-Cai-DasGupta 2001). Wilson for single proportions; Newcombe-Wilson or MN for differences/ratios.
- **Miettinen-Nurminen score CI is the regulatory standard for RD/RR** -- consistent with Pearson chi-square; preferred for NI margins.
- **Bender 2002 NNT convention:** when RD CI crosses zero, NNT CI is disjoint (NNTB-infinity-NNTH). Report as "NNTB X (NNTB Y to inf to NNTH Z)" -- standard in BMJ/Lancet/Cochrane.
- **NNT requires baseline risk to be meaningful** -- same OR gives dramatically different NNTs at different baseline risks.
- **OR overstates RR when prevalence > 10%.** Use modified Poisson (Zou 2004) with HC1/HC3 sandwich SE to estimate RR directly. Log-binomial is an alternative but often fails to converge.
- **HC3 sandwich is recommended for n <=250** (Long-Ervin 2000). HC1 (Stata default) and HC3 (R sandwich default) can differ enough to flip NI p-values.
- **Hauck-Donner effect** (small cell, large coefficient, small Wald chi-square) requires profile-likelihood inference. Detect with R `VGAM::hdeff()`.
- **Cornfield 1956 exact OR CI is over-conservative.** Mid-p (Berry-Armitage 1995) is less conservative; profile likelihood is transformation-invariant.
- **Forest plot reciprocal effects (OR 0.5 vs 2.0) are equidistant on log scale only.** Always use log scale.

## Related Skills

- clinical-biostatistics/logistic-regression - Adjusted OR + g-computation (in depth)
- clinical-biostatistics/categorical-tests - 2x2 testing producing the ORs
- clinical-biostatistics/subgroup-analysis - Forest plots and stratified effects
- clinical-biostatistics/survival-analysis - HR and RMST as effect measures
- clinical-biostatistics/trial-reporting - CONSORT 2025 + ICH E9(R1) effect reporting
- clinical-biostatistics/missing-data-sensitivity - Effect-measure CIs under MI pooling
- clinical-biostatistics/multiplicity-graphical - Effect reporting post-multiplicity adjustment
- machine-learning/survival-analysis - Predictive survival
