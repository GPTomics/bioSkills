# Trial Reporting - Usage Guide

## Overview

Prepares statistical reports for clinical trials following CONSORT 2025 (Hopewell et al *Lancet* 2025), SPIRIT 2025, ICH E9(R1) Addendum on Estimands (2019), and FDA 2023 Final Guidance "Adjusting for Covariates in RCTs." Covers Table 1 generation, analysis populations (ITT vs FAS vs PP vs Safety), the 5 ICH E9(R1) intercurrent-event strategies, MMRM under MAR with Kenward-Roger via mmrm package, reference-based MI via rbmi (J2R/CR/CIR), Permutt tipping-point sensitivity, and the Cro vs Bartlett information-anchored vs frequentist variance debate.

## Prerequisites

```bash
pip install tableone statsmodels scikit-learn pandas numpy
```

R is **strongly recommended** for confirmatory regulatory work (Python lacks Kenward-Roger):

```r
install.packages(c('mmrm', 'rbmi', 'gMCP', 'rpact', 'gsDesign'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Define the ICH E9(R1) estimand for my primary endpoint (5 attributes including ICE strategy)"
- "Fit MMRM with Kenward-Roger via the mmrm R package; UN covariance with pre-specified fallback hierarchy"
- "Run J2R reference-based MI via rbmi as MNAR sensitivity for treatment-policy estimand; report BOTH Rubin's and CMI+jackknife variance"
- "Compute Permutt tipping-point delta that flips primary p > 0.05; report in residual SD units"
- "Generate Table 1 with TableOne; pre-specify CONSORT 2025 + SPIRIT 2025 reporting per item 4 (data sharing) + item 21c (missing data) + item 15 (harms)"

## Example Prompts

### Estimand specification

> "My primary endpoint is change in HbA1c at 24 weeks. ICE is treatment discontinuation due to AE. Define the ICH E9(R1) estimand with 5 attributes and choose between treatment policy, hypothetical, composite, while-on-treatment, or principal stratum."

> "Cite Olarte Parra 2022 to argue MMRM-MAR IS a causal hypothetical estimator under specific identifying assumptions. Articulate the hypothetical scenario in the SAP."

### MMRM with Kenward-Roger

> "Implement MMRM in R mmrm package with UN covariance + Kenward-Roger-Linear (matches SAS PROC MIXED bit-for-bit). Pre-specify fallback: UN+KR → UN+Satterthwaite → heterogeneous Toeplitz → AR(1) → CS."

### Reference-based MI

> "Apply J2R imputation in rbmi for treatment-discontinuation ICEs in my treatment-policy estimand. Use both Bayesian MI + Rubin's rules (Cro 2019 information-anchored) AND CMI+jackknife (Wolbers 2022 frequentist). Justify discrepancy in CSR."

### Tipping-point sensitivity

> "Permutt 2016 tipping-point analysis: apply delta only to active-arm imputed values; scan delta from 0 to 20 in units of residual SD. Report the tipping delta that flips primary p > 0.05. Judge clinical plausibility against MCID."

### Aducanumab/Aprocitentan-style differential dropout

> "My trial has 35% dropout in active arm vs 12% in placebo, predominantly due to AEs. MAR is implausible. Switch from MMRM-MAR primary to treatment-policy with retrieved-dropout MI as primary, J2R as sensitivity (Aprocitentan 2024 precedent)."

### Table 1 and analysis populations

> "Generate Table 1 with age, sex, race, BMI, baseline severity. Group by ARM. Use SMD (NOT p-values) to assess balance. Distinguish ITT, FAS, PP, Safety populations with explicit pre-specified criteria."

### CONSORT 2025 reporting

> "Review my analysis against CONSORT 2025 (Hopewell 2025 *Lancet* 405:1633). Verify items 4 (data/code sharing), 8 (PPI), 15 (harms), 21c (missing data), 24a/24b (TIDieR intervention description)."

> "Estimands did NOT make consensus for mandatory CONSORT 2025 inclusion (Box 1 terminology only). Report the estimand per ICH E9(R1) directly."

## What the Agent Will Do

1. Define the ICH E9(R1) estimand with 5 attributes BEFORE choosing analysis method
2. Examine DS (Disposition) domain to characterise dropout patterns by arm
3. Choose estimand strategy based on clinical reasoning (treatment policy / hypothetical / composite / WoT / principal stratum)
4. Generate Table 1 with SMD-based balance assessment
5. Define ITT, FAS, PP, Safety populations with explicit pre-specified criteria
6. Execute MMRM-MAR (R mmrm) OR reference-based MI (R rbmi) per pre-specified estimand
7. Run MNAR sensitivity analyses including Permutt tipping-point in residual SD units
8. Apply graphical multiplicity (Bretz-Maurer via gMCP) for primary + key secondary
9. Format reporting per CONSORT 2025 30-item checklist + ICH E9(R1) estimand statement

## Tips

- **The estimand comes FIRST.** Kahan 2023 *Am J Epidemiol* 192:987: 98% of published trials don't articulate the estimand. Define 5 ICH E9(R1) attributes BEFORE choosing the method.
- **Aducanumab (2021)** is the textbook case where MAR-primary in trial with high differential missingness was regulator-divisive (6-1 AdCom against). Switch to treatment-policy + reference-based MI when differential dropout suggests informative missingness.
- **Aprocitentan (2024) precedent** is now the de facto FDA standard: hybrid imputation -- J2R for treatment-discontinuation ICEs, MMRM-MAR for other missingness.
- **Wegovy/Ozempic STEP precedent:** retrieved-dropout MI for treatment-policy when post-ICE data collected (FDA 2025 obesity guidance endorses).
- **MMRM with Kenward-Roger** is the FDA-favoured continuous-endpoint analysis. `mmrm` package `method = "Kenward-Roger-Linear"` matches SAS PROC MIXED bit-for-bit. Pre-specify the convergence fallback hierarchy in the SAP.
- **Reference-based MI variance debate (Cro vs Bartlett) is unsettled.** Report both Rubin's (information-anchored per Cro 2019) AND frequentist (CMI+jackknife per Wolbers 2022) for safety.
- **Permutt tipping-point delta should be in residual SD units** (FDA preference) for cross-trial comparison.
- **LOCF is biased even under MCAR** (Mallinckrodt 2008). NRC 2010 Rec 11 rejects LOCF. Don't use as "conservative" sensitivity.
- **Selection models (Diggle-Kenward) should be sensitivity only.** FDA prefers pattern-mixture (reference-based MI) because it is clinically articulable.
- **`sklearn.IterativeImputer` caveats:** `sample_posterior=True` only works with BayesianRidge (silently ignored otherwise). For confirmatory work, use R `rbmi` or `mice`.
- **Imputation model must be congenial** with analysis model (Meng 1994). If analysis includes interactions, imputation must too.
- **FAS vs ITT distinction matters:** FAS may exclude eligibility failures or no-post-baseline. ITT cannot. Pre-specify both with explicit FAS exclusion criteria.
- **Senn 1994 baseline balance testing is incoherent.** Don't condition adjustment on observed imbalance -- destroys nominal Type-I error.
- **CONSORT 2025 (Hopewell 2025) is 30 items**, not 25. 7 new items, 3 revised. Estimands in Box 1 only (consensus failed for mandatory). No DOORS framework.
- **SPIRIT 2025 has 34 items** (was 33). New PPI item (11) and trial-monitoring item (29).

## Related Skills

- clinical-biostatistics/cdisc-data-handling - SDTM/ADaM data preparation
- clinical-biostatistics/logistic-regression - Primary analysis with marginal-vs-conditional
- clinical-biostatistics/effect-measures - Effect-measure CIs under MI pooling
- clinical-biostatistics/subgroup-analysis - Pre-specified subgroup analyses for CSR
- clinical-biostatistics/missing-data-sensitivity - MMRM/rbmi/tipping-point in depth
- clinical-biostatistics/multiplicity-graphical - Bretz-Maurer graphs for co-primary
- clinical-biostatistics/survival-analysis - Time-to-event estimand framing
- clinical-biostatistics/power-and-sample-size - Sample size justification per CONSORT 2025 item 7
- clinical-biostatistics/adaptive-designs - Adaptive trial reporting
- clinical-biostatistics/bayesian-trials - Bayesian inference reporting
- reporting/rmarkdown-reports - Formatted statistical report generation
- experimental-design/multiple-testing - General multiplicity correction
