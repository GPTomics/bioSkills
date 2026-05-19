# Bayesian Trials - Usage Guide

## Overview

Designs Bayesian clinical trials including Phase I dose-finding (BOIN, CRM, EWOC, mTPI-2), meta-analytic-predictive (MAP) priors with robust mixtures for external data borrowing, EXNEX for basket trials, hierarchical models for safety AE multiplicity (Berry-Berry), Bayesian platform trials (I-SPY 2, GBM AGILE, REMAP-CAP), and posterior probability stopping rules. Covers FDA Bayesian Devices Guidance (2010), FDA Bayesian Methodology in Drugs Draft (January 2026), BOIN Fit-for-Purpose qualification (December 2021), and Project Optimus dose-optimisation.

## Prerequisites

R is the regulatory de facto standard:

```r
install.packages(c('RBesT', 'OncoBayes2', 'BOIN', 'dfcrm', 'escalation',
                   'trialr', 'bayesDP', 'psborrow2', 'rstan', 'cmdstanr',
                   'brms', 'bhmbasket', 'c212'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Phase 1 oncology BOIN design for 6 dose levels, target DLT 30%, max 30 patients"
- "Build a robust MAP prior in RBesT from 4 historical control arms for my Phase 3"
- "EXNEX basket trial across 5 tumour-type strata with 0.5 EX / 0.5 NEX mixture"
- "Pediatric extrapolation via power prior with discount gamma=0.5 per FDA Jan 2026 draft"
- "Posterior probability stopping rule for adaptive Bayesian trial; calibrate threshold via simulation to control frequentist Type-I at 0.025"

## Example Prompts

### Phase 1 dose-finding

> "Design BOIN for Phase 1 oncology: 6 dose levels, target DLT 0.30, cohort size 3, max 30 patients. Generate the escalation table for the protocol and simulate operating characteristics over true DLT rates (0.05, 0.10, 0.20, 0.30, 0.40, 0.55)."

> "CRM with Lee-Cheung 2009 calibrated skeleton. Target 0.30, 5 doses, n=30, starting dose 1. Compare OCs to BOIN."

### MAP priors

> "Fit a MAP prior in RBesT to 4 historical control arms (binomial outcome). Add robust mixture component with weight 0.2 for prior-data conflict protection. Report effective sample size of the prior."

> "Compare gMAP() with halfnormal tau prior 0.25 vs 0.50. Show effect on ESS and posterior."

### EXNEX basket trial

> "Design an EXNEX basket trial for 5 tumour-type strata. Default 0.5 EX / 0.5 NEX. Sensitivity over mixture weights (0.1, 0.3, 0.5, 0.7, 0.9). Implement in OncoBayes2."

### Pediatric extrapolation

> "Apply a power prior with gamma=0.5 to borrow 50% of adult trial data weight for my pediatric extrapolation. Implement in psborrow2 (FDA-supported)."

### Adaptive Bayesian stopping

> "Adaptive trial with 3 interim looks. Stop for efficacy if posterior P(theta > 0 | data) > 0.975. Calibrate this threshold via simulation under null to control frequentist Type-I at 0.025 one-sided."

### I-SPY 2 platform

> "Implement I-SPY 2 graduation criterion: arm graduates to Phase 3 if posterior predictive probability of success in 300-patient Phase 3 trial >= 0.85."

### Safety AE hierarchical

> "Berry-Berry 3-level hierarchical model for AE multiplicity across MedDRA PTs within SOCs. Spike-and-slab prior on log OR. Run in JMP Clinical or via c212 R package."

## What the Agent Will Do

1. Identify the Bayesian context (dose-finding, borrowing, basket, platform, sensitivity, safety)
2. Choose the appropriate model (BOIN, CRM, MAP, EXNEX, hierarchical)
3. Specify prior with sensitivity over hyperparameters
4. Simulate frequentist operating characteristics (Type-I, power) under regulatory expectation
5. Generate convergence diagnostics (R-hat <1.01, ESS >1000 per chain)
6. Validate against pinned package versions; document Docker/renv environment for reproducibility

## Tips

- **FDA preferred Bayesian device guidance is 2010** (CDRH). The first DRUG-side Bayesian guidance is the January 2026 CDER draft (FDA-2025-D-3217); comment period closed March 2026.
- **BOIN is FDA Fit-for-Purpose qualified (December 2021)** -- first dose-finding design with formal FDA endorsement. Investigator uses pre-tabulated escalation table; no bedside Bayesian software.
- **BOIN vs CRM:** BOIN is operationally simpler and more transparent; CRM is statistically more efficient under correct skeleton. The choice depends on team Bayesian sophistication and regulatory familiarity.
- **MAP priors require sensitivity over tau prior and mixture weight.** RBesT's robustify() with weight 0.2 is standard; verify with posterior predictive checks for prior-data conflict.
- **EXNEX (Neuenschwander 2016)** avoids hierarchical-model "catastrophic borrowing" when one basket truly differs. Default 0.5/0.5 weights with sensitivity analysis is standard.
- **Hemmings-Koch 2019 critique of Bayesian shrinkage:** appropriate for *replication planning*, NOT for *signal generation*. Bayesian shrinkage pre-emptively damps the heterogeneity you're searching for.
- **Power priors with gamma in 0.3-0.6** is the FDA Jan 2026 draft default for pediatric extrapolation. gamma = 1 (full pooling) ignores between-study heterogeneity.
- **Posterior probability stopping does not control frequentist Type-I automatically.** Calibrate threshold via simulation under null to demonstrate frequentist properties (ICH E20 expectation).
- **I-SPY 2 graduation criterion PP >= 0.85** is the canonical Bayesian platform standard (Barker 2009).
- **Spiegelhalter 1986 skeptical/enthusiastic prior framework** is still cited in modern Bayesian trial protocols for regulatory sensitivity.
- **Reproducibility for FDA submissions:** Docker container + renv-pinned R + Stan version; seeds + R-hat <1.01 + ESS >1000-2000 per parameter.
- **Project Optimus + FDA Bayesian Jan 2026 draft together signal that model-assisted Bayesian designs are the default for early oncology** as of 2026.

## Related Skills

- clinical-biostatistics/adaptive-designs - Group-sequential, SSR, platform trials
- clinical-biostatistics/subgroup-analysis - Bayesian shrinkage (Dixon-Simon, Berry)
- clinical-biostatistics/power-and-sample-size - Predictive probability of success
- clinical-biostatistics/multiplicity-graphical - Berry-Berry AE hierarchical
- clinical-biostatistics/trial-reporting - Bayesian inference reporting
- clinical-biostatistics/missing-data-sensitivity - Bayesian rbmi
- machine-learning/biomarker-discovery - Bayesian HTE
- experimental-design/sample-size - General methods
