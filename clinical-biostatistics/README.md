# clinical-biostatistics

## Overview

Postdoc/regulatory-grade statistical analysis methods for clinical trial data, spanning CDISC SDTM/ADaM data preparation, ICH E9(R1) estimand framework, FDA 2023 covariate adjustment, modern multiplicity (Bretz-Maurer graphical), survival analysis under non-PH (RMST, Fine-Gray vs cause-specific Cox), missing-data sensitivity (MMRM, reference-based MI, Permutt tipping point), adaptive designs (ICH E20 draft, BOIN, Mehta-Pocock), and Bayesian platform trials (I-SPY 2, EXNEX, MAP priors). Compliant with CONSORT 2025, SPIRIT 2025, FDA Multiple Endpoints 2022 Final, FDA NI 2016, and FDA Bayesian Methodology January 2026 draft.

**Tool type:** mixed | **Primary tools:** statsmodels, scipy, tableone, pyreadstat; R packages mmrm, rbmi, gMCP, rpact, gsDesign, RBesT, BOIN, survival, survRM2, riskRegression, lifelines

## Skills

| Skill | Description |
|-------|-------------|
| cdisc-data-handling | Read CDISC SDTM/ADaM (ADSL/BDS/OCCDS/ADTTE) with traceability, Pinnacle 21 validation, Define-XML 2.1, Dataset-JSON (Dec 2025), and the ADaMIG v3.0 roadmap |
| logistic-regression | Binary/ordinal/multinomial models with FDA 2023 marginal-vs-conditional estimand, g-computation, Brant test, Firth penalty, Hauck-Donner detection, modified Poisson for RR |
| categorical-tests | Chi-square, Boschloo's exact (uniformly more powerful than Fisher), CMH, mid-p McNemar (Fagerland 2013), Wilson/Newcombe/Miettinen-Nurminen CIs |
| effect-measures | OR/RR/RD/HR/NNT with calibrated CIs (Newcombe-MN-MOVER-profile-likelihood), Bender NNT convention, marginal vs conditional non-collapsibility, FDA 2023 g-computation |
| subgroup-analysis | Mantel-Haenszel, interaction tests, RERI, modern HTE (STEPP, SIDES, causal forests, X/R-learners), Bayesian shrinkage (Dixon-Simon, EXNEX), EMA 2019 subgroup guideline, gMCP, Yadlowsky RATE 2025 |
| trial-reporting | Table 1, ITT/FAS/PP/Safety, ICH E9(R1) 5 estimand strategies, Permutt tipping point, MMRM under MAR, reference-based MI (J2R/CR/CIR per Carpenter-Roger 2013), CONSORT 2025 + SPIRIT 2025 |
| survival-analysis | Cox PH + Therneau-Grambsch diagnostics, RMST (Royston-Parmar/Uno), Fine-Gray vs cause-specific Cox (Andersen-Keiding critique), MaxCombo + Magirr-Burman directional constraints, recurrent events (AG/PWP/WLW), interval censoring |
| missing-data-sensitivity | MMRM with Kenward-Roger (mmrm), reference-based MI via rbmi (Wolbers 2022), Cro vs Bartlett information-anchored vs frequentist variance debate, Permutt 2016 tipping-point, NRC 2010 framework |
| power-and-sample-size | Schoenfeld 1981 + Lakatos 1988 for non-PH, Fleiss continuity-correction debate, FDA 2016 NI double discount, TOST/Schuirmann, MCID vs δ distinction (Jaeschke 1989, Norman 0.5 SD heuristic) |
| multiplicity-graphical | Bretz-Maurer graphical procedures via gMCP, gatekeeping (parallel/serial/mixed), closed-testing principle (Marcus-Peritz-Gabriel; Goeman 2021 admissibility), FDA Multiple Endpoints October 2022 Final |
| adaptive-designs | Group-sequential (rpact, gsDesign), blinded/unblinded SSR, Mehta-Pocock promising zone, RAR ethics (Hey-Kimmelman vs Berry; Robertson 2023 review), Bayesian platforms (I-SPY 2, REMAP-CAP), ICH E20 Step 2b/3 draft (June 2025) |
| bayesian-trials | BOIN (FDA Fit-for-Purpose Dec 2021), CRM/EWOC/mTPI-2, MAP + robust MAP via RBesT, EXNEX basket trials, Berry-Berry AE hierarchical, Spiegelhalter skeptical priors, FDA Bayesian Methodology Jan 2026 draft |

## Example Prompts

- "Load my CDISC .xpt files and create a subject-level analysis dataset; validate against Pinnacle 21"
- "Run logistic regression with marginal RD via g-computation per FDA 2023 covariate adjustment guidance"
- "Test treatment-outcome association with Boschloo's exact test (more powerful than Fisher) and report Miettinen-Nurminen CI"
- "Compute marginal vs conditional ORs for my RCT and explain non-collapsibility per Permutt 2020"
- "Apply causal forests for HTE detection with RATE/AUTOC omnibus test; report bias-corrected subgroup estimates"
- "Set up MMRM with mmrm and Kenward-Roger correction; add J2R sensitivity via rbmi with both Cro information-anchored and Wolbers frequentist variance"
- "Compute tipping-point delta in residual SD units that would flip my primary p-value above 0.05"
- "Fit Cox PH for OS; if PH fails, switch to RMST at pre-specified tau=36 months"
- "Design a graphical multiplicity procedure in gMCP for primary + 2 key secondary + 1 subgroup hypothesis"
- "Schoenfeld events formula for HR=0.70 OS; switch to Lakatos simulation under expected delayed-effect immunotherapy"
- "FDA 2016 NI margin: apply double discount with historical M1 lower-CI 20% to get M2=10%; verify M2 < MCID"
- "BOIN Phase 1 escalation table for 6 doses with target DLT 30%; FDA Fit-for-Purpose qualified design"
- "Robust MAP prior in RBesT from 4 historical control arms with weight 0.2 vague mixture"

## Requirements

```bash
pip install statsmodels scipy pingouin tableone pyreadstat pandas numpy matplotlib scikit-learn lifelines scikit-survival
```

Optional for Firth penalty:

```bash
pip install firthmodels
```

R is **strongly recommended** for confirmatory regulatory work:

```r
install.packages(c(
    'mmrm', 'rbmi',                                   # MMRM + reference-based MI
    'gMCP', 'graphicalMCP',                           # multiplicity
    'rpact', 'gsDesign', 'gsDesign2',                 # group-sequential + adaptive
    'survival', 'survRM2', 'cmprsk', 'riskRegression', 'mstate', 'flexsurv', 'icenReg',  # survival
    'BOIN', 'dfcrm', 'escalation', 'trialr',          # phase 1 dose-finding
    'RBesT', 'OncoBayes2', 'bayesDP', 'psborrow2',    # Bayesian
    'pwr', 'presize', 'npsurvSS',                     # power/sample size
    'mice', 'mitools'                                 # multiple imputation
))
```

## Regulatory Reference Documents

- **CONSORT 2025** (Lancet, April 2025) -- 30-item checklist for trial reporting
- **SPIRIT 2025** (Lancet, April 2025) -- 34-item protocol guideline
- **ICH E9(R1) Addendum on Estimands** (final November 2019; EMA effective Feb 2020; FDA May 2021)
- **ICH E20 Adaptive Clinical Trials** -- Step 2b/3 draft June 2025 (NOT final as of May 2026)
- **FDA Adjusting for Covariates in RCTs** (Final May 2023)
- **FDA Multiple Endpoints in Clinical Trials** (Final October 2022)
- **FDA Master Protocols** (Final March 2022)
- **FDA Adaptive Designs for Drugs and Biologics** (Final December 2019)
- **FDA Non-Inferiority Trials** (Final November 2016)
- **FDA Use of Bayesian Methodology in Clinical Trials** (Draft January 2026)
- **FDA Project Optimus** + Oncology Dose-Optimisation (Final August 2024)
- **FDA BOIN Fit-for-Purpose Qualification** (December 2021)
- **EMA Investigation of Subgroups in Confirmatory Trials** (effective August 2019)
- **EMA Adjustment for Baseline Covariates** (effective September 2015)
- **EMA Missing Data in Confirmatory Clinical Trials** (2010)
- **NRC 2010** *Prevention and Treatment of Missing Data in Clinical Trials* (National Academies)

## Related Skills

- **machine-learning** - Survival prediction, biomarker discovery, HTE estimation
- **experimental-design** - Power analysis, sample size, multiple testing (general methods)
- **epidemiological-genomics** - Genomic epidemiology and AMR surveillance
- **workflows** - Clinical trial analysis pipeline
- **reporting** - Quarto/RMarkdown report generation for CSR sections
