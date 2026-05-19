# Adaptive Clinical Trial Designs - Usage Guide

## Overview

Designs adaptive clinical trials with one or more pre-specified interim adaptations: group-sequential boundaries (O'Brien-Fleming, Pocock, Lan-DeMets spending), sample-size re-estimation (blinded Friede-Kieser, unblinded Cui-Hung-Wang, Mehta-Pocock promising zone), seamless Phase 2/3 with treatment-arm selection, population enrichment, and response-adaptive randomisation. Covers FDA 2019 Final Adaptive Designs Guidance, FDA 2022 Master Protocols, and ICH E20 Step 2b/3 draft (June 2025, NOT final).

## Prerequisites

R is the regulatory de facto standard:

```r
install.packages(c('rpact', 'gsDesign', 'gsDesign2', 'adaptr', 'simtrial',
                   'BOIN', 'dfcrm', 'escalation', 'trialr'))
```

Commercial software for confirmatory submissions: East/EastHorizon (Cytel), ADDPLAN (ICON), FACTS (Berry Consultants).

## Quick Start

Tell your AI agent what you want to do:
- "Design a group-sequential trial with O'Brien-Fleming boundaries in gsDesign for OS, HR=0.70, 90% power"
- "Implement blinded sample-size re-estimation for my continuous endpoint trial with uncertain variance"
- "Set up Mehta-Pocock promising-zone SSR with CHW weighting"
- "BOIN phase 1 dose-finding for 6 dose levels, target DLT 30%, cohort size 3"
- "Bayesian platform trial design with 4 arms, response-adaptive randomisation, futility stopping"

## Example Prompts

### Group-sequential

> "Pivotal trial OS endpoint. Design O'Brien-Fleming boundary in gsDesign with 3 interim looks at 33%/67%/100% information. HR=0.65, control hazard 0.04/month, 24-month accrual, 12-month minimum follow-up. 90% power at one-sided alpha 0.025."

> "Use Lan-DeMets sfLDOF spending function for flexibility in actual look timing."

### Blinded SSR

> "Continuous primary endpoint; expected SD 12 but uncertain. Implement blinded SSR at interim n=200 with sample-size recalculation up to n_max=300. No Type-I inflation per Friede-Kieser 2006."

### Mehta-Pocock promising zone

> "Pre-specify Mehta-Pocock promising-zone SSR. Initial n=300/arm; if interim conditional power ∈ (0.3, 0.8), increase to n_max=450/arm with CHW weights from rpact getDesignInverseNormal."

### Adaptive enrichment

> "Phase 3 trial with biomarker subgroup. Implement adaptive enrichment via closed testing: at interim, if biomarker-positive subgroup has CP > 0.7 and full population has CP < 0.3, drop the negative subgroup and re-power on positive."

### Bayesian platform (I-SPY 2 style)

> "Multi-arm biomarker-stratified Bayesian platform for breast cancer neoadjuvant. Implement RAR with graduation criterion: PP(success in 300-patient Phase 3) >= 0.85. Use FACTS or custom Stan."

### Phase 1 dose-finding

> "BOIN for 6 dose levels in Phase 1 oncology. Target DLT rate 30%, cohort size 3, max sample size 30. Generate the pre-tabulated escalation table for the protocol."

> "Compare CRM (1-parameter logistic skeleton) vs BOIN OCs for the same setting. Report MTD selection accuracy."

### Project Optimus

> "Project Optimus-compliant Phase 1b/2: after BOIN finds MTD, randomise to two doses (MTD and MTD/2). Multi-arm BOIN-12 design with efficacy + toxicity."

## What the Agent Will Do

1. Identify the adaptation type (early stopping, SSR, enrichment, treatment selection, RAR)
2. Pre-specify the design in rpact or gsDesign
3. Compute operating characteristics via simulation (Type-I, power, expected sample size)
4. Choose the boundary (OBF, Pocock, Wang-Tsiatis) or spending function (Lan-DeMets)
5. For SSR: specify CHW weights or promising-zone CP boundaries
6. For RAR: pre-specify time-trend covariates and analysis weights
7. Generate the SAP-ready specification (design object, operating characteristics, simulation script)

## Tips

- **FDA prefers group-sequential over RAR** for confirmatory trials. RAR carries drift bias, estimator bias, and operational complexity that FDA finds unacceptable in 2-arm confirmatory settings.
- **OBF is conservative early, near-nominal final** (alpha ~0.024 of 0.025 at final analysis with k=4) -- preferred by FDA over Pocock (which spends alpha equally and has larger final penalty).
- **Lan-DeMets spending function** allows analyses at different information fractions than originally planned -- operational flexibility valued by FDA.
- **Blinded SSR (Friede-Kieser 2006) is uncontroversial** because Type-I is not affected. Unblinded SSR (Cui-Hung-Wang 1999) requires CHW weighting to preserve Type-I.
- **Mehta-Pocock promising zone is calibrated** so Type-I inflation is negligible (~0.001) without CHW weighting -- but the Jennison-Turnbull 2015 critique says it's "stealth alpha inflation." Pre-specify the increase rule transparently and report simulation operating characteristics.
- **DMC firewall is essential** for unblinded SSR. Sponsor receives only the increase/no-increase decision; interim effect estimate never leaks. ICH E20 codifies this.
- **ICH E20 is NOT final** as of May 2026. It is Step 2b/3 draft (June 25 2025). Step 4 targeted late 2026. Do not cite as final guidance.
- **FDA Master Protocols Guidance was finalised March 2022**, NOT 2018 (2018 was the draft).
- **BOIN was FDA Fit-for-Purpose qualified December 2021** under the Drug Development Tools program -- first formal FDA endorsement of a specific dose-finding design.
- **Project Optimus (FDA OCE, 2021; final dose-optimisation guidance Aug 2024)** rewrote Phase 1 oncology -- randomised dose comparison before registration, replacing MTD-and-go.
- **RAR consensus (Robertson 2023 Stat Sci 38:185):** appropriate for multi-arm (≥3 arms), rare disease, biomarker-stratified; INAPPROPRIATE for 2-arm confirmatory.
- **For RAR analysis, pre-specify time-trend covariates** -- drift bias is the largest operational risk.
- **Bayesian primary inference in pivotals is permitted by FDA Bayesian Jan 2026 draft** provided simulation-based Type-I error calibration.

## Related Skills

- clinical-biostatistics/power-and-sample-size - Sample size for adaptive designs
- clinical-biostatistics/multiplicity-graphical - Closed testing in adaptive contexts
- clinical-biostatistics/bayesian-trials - Bayesian platform, BOIN/CRM/EWOC
- clinical-biostatistics/trial-reporting - Reporting per CONSORT 2025
- clinical-biostatistics/survival-analysis - Adaptive TTE
- experimental-design/sample-size - General methods
