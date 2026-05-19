# Missing Data Sensitivity - Usage Guide

## Overview

Implements pre-specified missing-data analyses for confirmatory clinical trials following NRC 2010 and ICH E9(R1). Covers MMRM under MAR with Kenward-Roger correction, reference-based multiple imputation (J2R, CR, CIR, LMCF per Carpenter-Roger 2013), Permutt delta-adjustment / tipping-point analysis, pattern-mixture identifying restrictions, and the Cro vs Bartlett information-anchored vs frequentist variance debate.

## Prerequisites

```bash
pip install scikit-learn statsmodels pandas numpy
```

R is strongly recommended for confirmatory work (the Python alternatives lack Kenward-Roger):

```r
install.packages(c('mmrm', 'rbmi', 'mice', 'mitools'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Fit an MMRM under MAR for continuous longitudinal endpoint with Kenward-Roger correction"
- "Run J2R reference-based multiple imputation as MNAR sensitivity for my treatment-policy estimand"
- "Compute the Permutt tipping-point delta that flips the primary p-value above 0.05"
- "Report both Cro information-anchored variance AND Wolbers frequentist variance for J2R analysis"
- "Convert my CDISC trial dataset to long format and apply MMRM with the SAS PROC MIXED-equivalent options"

## Example Prompts

### MMRM under MAR

> "Fit an MMRM in R `mmrm` for change-from-baseline outcome with treatment, visit, treatment-by-visit, baseline, baseline-by-visit. Use unstructured covariance with Kenward-Roger-Linear DF. Report the contrast at the primary timepoint."

> "MMRM with UN+KR failed to converge. Fall back to heterogeneous Toeplitz with KR, per pre-specified SAP hierarchy."

### Reference-based MI

> "My estimand is treatment policy. ICE = treatment discontinuation due to AE in active arm. Run J2R imputation in rbmi with Bayesian MI + Rubin's rules as primary, CMI+jackknife frequentist as supportive."

> "Apply CIR (copy-increments-in-reference) imputation rather than J2R because we want to credit the on-treatment increment but assume post-ICE trend matches placebo."

### Tipping point

> "Run a tipping-point analysis. Apply delta only to active-arm imputed values, scan delta from 0 to 20 in units of the residual SD. Report the tipping delta that flips primary p above 0.05."

> "Tipping delta = 8 mmHg. Compare to MCID of 5 mmHg for SBP -- judge clinical plausibility."

### Variance debate

> "Report both Rubin's rules variance (Cro 2019 information-anchored) and CMI+jackknife frequentist variance (Wolbers 2022) for the J2R primary sensitivity analysis. Justify the discrepancy in the report."

### Aducanumab/Aprocitentan-style scenarios

> "My trial has differential dropout -- 35% in active arm vs 12% in placebo, with active-arm dropouts predominantly due to AEs. MAR is implausible. Switch to treatment-policy estimand with retrieved-dropout MI as primary, J2R as sensitivity."

## What the Agent Will Do

1. Examine DS (Disposition) domain or equivalent to characterise dropout patterns by arm
2. Determine if MAR is plausible (symmetric dropout) or implausible (differential)
3. Pre-specify the estimand strategy per ICH E9(R1) based on clinical reasoning
4. Execute the primary analysis under the chosen assumption (MMRM-MAR or reference-based MI)
5. Run pre-specified MNAR sensitivity analyses (J2R, CR, CIR, tipping-point)
6. Report both information-anchored (Rubin's) and frequentist (jackknife) variance for reference-based MI
7. Generate tipping-point report with delta in residual SD units (FDA preference)

## Tips

- **MAR vs MNAR is fundamentally untestable from observed data.** Pre-specify based on clinical reasoning (DS domain examination), not observed data.
- **LOCF is biased even under MCAR** (Mallinckrodt 2008). NRC 2010 Rec 11 explicitly rejects LOCF as default. Never use as primary or as "conservative" sensitivity.
- **MMRM under MAR is the FDA-favoured continuous-endpoint analysis** with UN covariance + Kenward-Roger. Use `mmrm` package in R; statsmodels.mixedlm lacks KR and is not FDA-equivalent.
- **Always pre-specify the convergence fallback hierarchy** in the SAP (UN+KR → UN+Satterthwaite → het Toeplitz → AR(1) → CS). Document any deviation invoked at analysis time.
- **Reference-based MI (Carpenter-Roger 2013) operationalises MNAR as clinical narrative**, not arbitrary delta. J2R, CR, CIR, LMCF are the four operational forms -- choose the one matching the clinical scenario.
- **The variance debate is ongoing.** EMA tolerates either Rubin's or frequentist; FDA increasingly wants frequentist supplement. **Report both** for safety.
- **Tipping-point delta should be in residual SD units** (FDA preference) for cross-trial comparison, not raw outcome units.
- **sklearn IterativeImputer has critical gotchas**: `sample_posterior=True` only works with BayesianRidge; silently ignored otherwise. For confirmatory work, use R `rbmi` or `mice`.
- **Imputation model must be congenial with analysis model** (Meng 1994). If analysis includes treatment-by-covariate interaction, imputation must too.
- **Selection models (Diggle-Kenward 1994)** should be sensitivity-only. FDA pushes back because MAR-vs-MNAR conclusion is driven by joint normality assumption, not data.

## Related Skills

- clinical-biostatistics/trial-reporting - Estimand framework + CONSORT 2025 reporting
- clinical-biostatistics/logistic-regression - Adjusted logistic with MI
- clinical-biostatistics/effect-measures - Effect-measure CIs under pooling
- clinical-biostatistics/cdisc-data-handling - DS domain reasoning for missingness
- clinical-biostatistics/survival-analysis - Informative censoring analogues
- clinical-biostatistics/adaptive-designs - Interim missing-data
