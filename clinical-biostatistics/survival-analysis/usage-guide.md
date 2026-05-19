# Survival Analysis - Usage Guide

## Overview

Performs time-to-event analysis for clinical trials with calibration to whether proportional hazards holds, whether competing events exist, whether censoring is informative, and which ICH E9(R1) estimand the trial targets. Covers Cox PH with Therneau-Grambsch diagnostics, restricted mean survival time (RMST) as a hazard-free alternative, Fine-Gray vs cause-specific Cox for competing risks, MaxCombo and weighted log-rank for non-proportional hazards, recurrent events (Andersen-Gill, PWP, WLW), and interval-censored data.

## Prerequisites

```bash
pip install lifelines scikit-survival statsmodels pandas numpy matplotlib
```

R is strongly recommended for production survival analysis:

```r
install.packages(c('survival', 'survRM2', 'cmprsk', 'riskRegression', 'mstate',
                   'flexsurv', 'icenReg', 'rpsftm'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Run a Cox proportional hazards regression for OS in my oncology trial; check the PH assumption"
- "Estimate the restricted mean survival time difference at 36 months as the primary endpoint because hazards cross"
- "I have a competing risk: cardiac death and non-cardiac death. Fit cause-specific Cox for both and report the cumulative incidence function"
- "Apply MaxCombo for a delayed-effect immunotherapy trial with pre-specified directional constraint"
- "Convert ADTTE CNSR=0/1/2 to the standard event indicator for lifelines"

## Example Prompts

### Standard Cox PH

> "Fit a Cox model for OS in my trial with treatment, age, and ECOG performance status. Run cox.zph and report the global PH test. Stratify on randomisation strata (region)."

> "Compute the hazard ratio for treatment vs placebo with 95% CI. If PH is violated, switch to RMST."

### RMST under non-PH

> "PH is clearly violated in my immunotherapy trial. Compute the RMST difference at 36 months as the primary endpoint with delta-method CI. Pre-specified tau is 36 from the SAP."

> "Compare the time-to-RMST tradeoff at tau=12, 24, 36, and 60 months for sensitivity."

### Competing risks

> "I have a cardiovascular outcome trial with cardiac death (event of interest) and non-cardiac death (competing). Fit cause-specific Cox for cardiac death AND for non-cardiac death. Report both cause-specific HRs and the cumulative incidence functions via Aalen-Johansen."

> "My oncology trial has progression and death. Compute the CIF for both via Aalen-Johansen and report PFS event-free survival treating either as event."

### Non-proportional hazards

> "Run MaxCombo with G(0,0), G(0,1), G(1,0), G(1,1) family for the late-effect immunotherapy. Pre-specified directional constraint: require positive z-statistic at the late-emphasis G(0,1) weight before claiming superiority."

> "Fit a Royston-Parmar flexible parametric model with 5-knot cubic spline for log cumulative hazard. Report HR(t) curves over follow-up."

### Recurrent events

> "Subjects have multiple hospitalisations during follow-up. Fit Andersen-Gill with robust SE clustered on USUBJID. Report rate ratio."

> "Hospitalisations may differ by event order (first vs recurrent). Fit PWP gap-time with stratum-specific baseline hazards."

### Interval censoring

> "PFS scans were every 12 weeks but the actual progression date is unknown. Fit interval-censored Cox via icenReg::ic_par to avoid the midpoint imputation bias."

### ADTTE handling

> "Load my ADTTE.csv with CNSR per CDISC convention (0=event, 1+=censoring reasons). Convert to event/time for lifelines and fit Cox."

## What the Agent Will Do

1. Load the time-to-event dataset (ADTTE preferred; convert CNSR if needed)
2. Validate KM curves are biologically plausible; tabulate censoring by arm
3. Fit Cox PH and run Therneau-Grambsch diagnostics (cox.zph in R, proportional_hazard_test in lifelines)
4. If PH holds: report HR with stratified log-rank; if PH violated: switch to RMST or time-varying Cox
5. For competing risks: report both cause-specific Cox AND Fine-Gray CIF; never use KM
6. For immuno-oncology with delayed effect: use MaxCombo with pre-specified directional constraint
7. Report effect estimates with 95% CI, p-values, and the ICH E9(R1) estimand the analysis targets

## Tips

- **cox.zph p > 0.05 does NOT prove PH.** Use it as a failure detector, not a validator. Always plot scaled Schoenfeld residuals graphically.
- **Pre-specify tau for RMST in the SAP.** Post-hoc tau tuning is p-hacking. Constrain tau <= min(largest follow-up across arms) to avoid extrapolation.
- **Never use KM in competing-risk settings.** It biases the survival estimate upward. Use Aalen-Johansen cumulative incidence.
- **Fine-Gray is for CIF prediction, NOT for causal inference.** Andersen-Keiding 2012 showed it violates the three principles for valid hazard functionals. Use cause-specific Cox for etiology.
- **MaxCombo can reject in opposite directions on the same data** (Magirr-Burman 2021 KEYNOTE-042 demonstration). Always pre-specify directional constraints.
- **ADTTE CNSR convention is the opposite of statistical packages.** CDISC: CNSR=0 means event; R/Python: event=1 means event. Always convert explicitly.
- **PFS is really interval-censored** but conventionally treated as right-censored at midpoint. This is acceptable when scan intervals are short and balanced; switch to interval-censored Cox when intervals > 4 weeks or differ between arms.
- **Stratification factors from randomisation MUST appear in analysis** (Kahan-Morris 2012). Use stratified log-rank and `strata()` in Cox, or include as covariates.
- **The Schoenfeld 1981 formula assumes PH.** Under non-PH (immuno-oncology), it under-estimates required events by 20-50%. Use Lakatos 1988 or simulation under the expected hazard pattern.
- **The Cox HR under PH violation is a time-averaged log-HR** (Xu-O'Quigley 2000). For a clinically interpretable summary, use RMST.

## Related Skills

- clinical-biostatistics/effect-measures - HR vs RMST as effect measures
- clinical-biostatistics/trial-reporting - ICH E9(R1) estimand framework for time-to-event
- clinical-biostatistics/missing-data-sensitivity - Informative censoring sensitivity
- clinical-biostatistics/cdisc-data-handling - ADTTE structure and CNSR convention
- clinical-biostatistics/subgroup-analysis - Subgroup HTE for survival endpoints
- clinical-biostatistics/power-and-sample-size - Schoenfeld and Lakatos for TTE
- clinical-biostatistics/multiplicity-graphical - Co-primary survival endpoints
- machine-learning/survival-analysis - Predictive survival models
