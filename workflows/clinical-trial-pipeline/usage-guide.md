# Clinical Trial Analysis Pipeline - Usage Guide

## Overview

End-to-end workflow for analyzing clinical trial data following CONSORT 2025, ICH E9(R1) estimand framework, and FDA 2023 covariate adjustment guidance. Takes CDISC SDTM/ADaM domain files through data preparation, ICH E9(R1) estimand specification, primary analysis (logistic regression with marginal vs conditional reporting, MMRM for continuous longitudinal, Cox/RMST for time-to-event), modern categorical testing, subgroup analysis with HTE methods, missing-data sensitivity (MMRM under MAR, reference-based MI, Permutt tipping point), graphical multiplicity (Bretz-Maurer), and regulatory-compliant reporting.

## Prerequisites

```bash
pip install statsmodels scipy tableone pyreadstat pandas numpy matplotlib scikit-learn lifelines scikit-survival
```

Optional for rare events:
```bash
pip install firthmodels
```

R is strongly recommended for confirmatory regulatory work (MMRM with Kenward-Roger; reference-based MI; graphical multiplicity):
```r
install.packages(c('mmrm', 'rbmi', 'gMCP', 'rpact', 'gsDesign', 'survival', 'survRM2', 'RBesT', 'BOIN'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Analyze my clinical trial data from start to finish with ICH E9(R1) estimand framework"
- "I have CDISC .xpt files for a vaccine trial -- run the full analysis with CONSORT 2025 reporting"
- "Perform a complete statistical analysis of treatment vs placebo with marginal RD primary per FDA 2023"
- "Run logistic regression, subgroup analysis with causal forests, and generate Table 1 for my trial data"
- "Continuous longitudinal endpoint: run MMRM with Kenward-Roger as primary, J2R reference-based MI as MNAR sensitivity"

## Example Prompts

### Full Pipeline

> "I have DM, AE, and EX domain files from a BCG vaccination trial. Define the ICH E9(R1) estimand first (5 attributes including ICE strategy). Create a subject-level dataset, run logistic regression on COVID-19 severity with marginal RD via g-computation per FDA 2023 (conditional OR as supportive), test associations with Boschloo's exact, analyze subgroups by patient load with graphical multiplicity via gMCP, and generate a Table 1 per CONSORT 2025."

> "Analyze my clinical trial data end to end. The primary endpoint is binary adverse event. I need: marginal RD via g-computation, conditional OR as supportive, subgroup forest plot with INTERACTION p-values (not per-subgroup comparison), Permutt tipping-point sensitivity analysis, and CONSORT 2025 flow diagram."

### Continuous longitudinal endpoint

> "Continuous primary endpoint is change in HbA1c at week 24. Pre-specify ICH E9(R1) hypothetical estimand for treatment discontinuation ICE. Fit MMRM with UN+Kenward-Roger via R mmrm. Sensitivity: J2R reference-based MI via R rbmi reporting BOTH Cro information-anchored variance AND Wolbers frequentist (CMI+jackknife)."

### Time-to-event primary endpoint

> "Primary endpoint is OS. Run stratified log-rank + Cox PH with Therneau-Grambsch diagnostic. If PH violated, switch to RMST at pre-specified tau=36 months. If competing risks (non-cardiac death), cause-specific Cox for both events + Aalen-Johansen CIF."

### Data Preparation

> "Load my CDISC SDTM files and create a subject-level analysis dataset. Merge demographics with adverse events, aggregate to one row per subject with maximum severity, code the treatment variable, and tabulate DS (Disposition) domain for differential dropout patterns by arm (informs missing-data strategy)."

### Primary Analysis

> "Run logistic regression on my prepared clinical dataset with treatment as the primary predictor, adjusting for age, sex, and randomisation strata. Extract BOTH conditional OR (model-fit, supportive) and marginal RD via g-computation with a percentile bootstrap CI per FDA 2023."

### Subgroup Analysis

> "Test whether the treatment effect varies across pre-specified subgroups via INTERACTION terms in single model (NOT per-subgroup p-comparisons). For continuous biomarker subgroups, use STEPP. Generate forest plot. Apply graphical multiplicity allocation via gMCP, allocating (for example) 20% of alpha to the subgroup family."

### Reporting

> "Generate Table 1 of baseline characteristics by arm with SMD (NOT p-values per Senn 1994 incoherence). Run missing-data sensitivity per CONSORT 2025 item 21c. Provide ICH E9(R1) estimand statement (5 attributes) in CSR appendix. Format harms per CONSORT 2025 item 15."

## What the Agent Will Do

1. Pre-specify ICH E9(R1) estimand (5 attributes including ICE strategy) BEFORE any analysis
2. Load CDISC SDTM/ADaM files (.xpt or .csv); convert ADTTE CNSR convention if needed
3. Tabulate DS (Disposition) domain dropout patterns by arm (informs missing-data strategy)
4. Aggregate event-level data to subject level; merge domains; verify USUBJID uniqueness
5. Generate Table 1 with SMD (NOT p-values)
6. Fit logistic regression with explicit reference category; report marginal RD via g-computation per FDA 2023 as primary
7. Run Boschloo (preferred over Fisher) for small 2x2 categorical tests; Pearson chi-square for adequate counts
8. Test subgroup interactions in single model; apply graphical multiplicity via gMCP
9. For continuous longitudinal: MMRM with Kenward-Roger as primary; reference-based MI as MNAR sensitivity
10. For time-to-event: Cox + cox.zph diagnostic; switch to RMST under PH violation
11. Permutt tipping-point sensitivity in residual SD units
12. Format reporting per CONSORT 2025 30-item checklist

## Tips

- Define the ICH E9(R1) estimand FIRST. Kahan 2023 Am J Epidemiol 192:987 documents 98% of trials don't articulate the estimand; pre-specify the 5 attributes in the SAP before choosing analysis.
- Per FDA 2023, primary estimand for binary endpoints is marginal RD via g-computation; conditional OR is a DIFFERENT parameter due to OR non-collapsibility (Permutt 2020).
- Boschloo's exact is uniformly more powerful than Fisher's exact at the same Type-I (Mehta-Senchaudhuri 2003); use as default for small 2x2.
- Stratified randomisation factors MUST appear in analysis (Kahan-Morris 2012); ignoring them makes inference conservative (SEs biased upward, Type-I below nominal, power lost).
- Always set explicit reference category in logistic regression (e.g., Placebo) to avoid alphabetical reversal of OR direction.
- Use INTERACTION terms in a single model to test subgroup effects; comparing per-subgroup p-values is statistically invalid.
- Report standardized mean differences (SMD > 0.1 as a notable-imbalance convention) rather than p-values for Table 1 balance; baseline hypothesis testing is incoherent (Senn 1994).
- Pre-specify subgroups before unblinding; post-hoc subgroup results are hypothesis-generating only per EMA 2019.
- When continuous longitudinal endpoint with monotone MAR: use R mmrm with method="Kenward-Roger-Linear" to match SAS PROC MIXED.
- For MNAR sensitivity with reference-based MI (J2R/CR/CIR per Carpenter-Roger 2013), report BOTH Rubin information-anchored variance AND frequentist CMI+jackknife (Cro vs Bartlett debate is unsettled).
- Permutt tipping-point delta should be in residual SD units (FDA preference) for cross-trial comparison.
- LOCF is biased even under MCAR (Mallinckrodt 2008; NRC 2010 Rec 11 rejects); never use as "conservative" sensitivity.
- For time-to-event, cox.zph p > 0.05 does NOT prove PH; use as failure detector with graphical residual plot as primary.
- Under PH violation, RMST is more interpretable than HR; pre-specify tau in SAP, NOT post-hoc.
- ADTTE CNSR convention is OPPOSITE of statistical packages (CDISC: CNSR=0 means event); always convert before passing to R/Python.
- For multiple endpoints, design a graphical procedure in R gMCP with pre-specified weights per FDA Multiple Endpoints Final Oct 2022.
- Run the pipeline in order. Data preparation issues (duplicate subjects, unmapped columns, DS not tabulated) cascade into every downstream step.

## Related Skills

- clinical-biostatistics/cdisc-data-handling - CDISC SDTM/ADaM, Pinnacle 21, ADTTE conventions
- clinical-biostatistics/logistic-regression - Marginal vs conditional, g-computation, Brant test, Firth
- clinical-biostatistics/categorical-tests - Boschloo, mid-p McNemar, Wilson/MN CIs
- clinical-biostatistics/effect-measures - NNT Bender 2002, profile likelihood, modified Poisson
- clinical-biostatistics/subgroup-analysis - Causal forests, STEPP, EXNEX, Yadlowsky RATE
- clinical-biostatistics/trial-reporting - ICH E9(R1) 5 estimand strategies, CONSORT 2025
- clinical-biostatistics/missing-data-sensitivity - MMRM/KR, J2R/CR/CIR, Permutt tipping point
- clinical-biostatistics/multiplicity-graphical - Bretz-Maurer graphs, closed-testing admissibility
- clinical-biostatistics/survival-analysis - Cox/RMST/Fine-Gray/MaxCombo/recurrent events
- clinical-biostatistics/power-and-sample-size - Schoenfeld/Lakatos, NI margin, MCID
- clinical-biostatistics/adaptive-designs - Group-sequential, SSR, platform trials
- clinical-biostatistics/bayesian-trials - BOIN, MAP priors, RWE
- reporting/rmarkdown-reports - Formatted statistical report generation
- reporting/publication-tables - Table 1 (SMD not baseline p-values, show missingness) and Word/LaTeX export
