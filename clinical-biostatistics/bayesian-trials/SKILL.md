---
name: bio-clinical-biostatistics-bayesian-trials
description: Designs Bayesian clinical trials including Phase I dose-finding (BOIN, CRM, EWOC, mTPI-2), meta-analytic-predictive (MAP) priors with robust mixtures for external data borrowing, EXNEX for basket trials, hierarchical models for safety AE (Berry-Berry), Bayesian platform trials (I-SPY 2, GBM AGILE, REMAP-CAP), and posterior probability stopping rules. Covers FDA Bayesian Devices Guidance (2010), FDA Bayesian Methodology in Drugs Draft (January 2026), BOIN Fit-for-Purpose qualification (December 2021), and Project Optimus dose-optimisation. Use when designing dose-finding studies, platform trials, or sensitivity analyses with informative priors.
tool_type: r
primary_tool: RBesT
goal_approach_exempt: true
---

## Version Compatibility

Reference examples tested with: R `RBesT` 1.7+ (Roche), `OncoBayes2` 0.8+ (Novartis), `BOIN` 2.7+, `dfcrm` 0.2-2+, `escalation` 0.1+, `trialr` 0.1.6+, `bayesDP`, `psborrow2` (FDA-supported), `rstan` / `cmdstanr`, `brms`. Legacy: `JAGS`, `WinBUGS`.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name`
- Confirmatory regulatory work: validate against pinned package versions in submission

If code throws an error, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Bayesian Clinical Trials

**"Design a Bayesian clinical trial"** -> Specify a prior, likelihood, and decision rule with frequentist operating characteristics demonstrated via simulation; for dose-finding use FDA-endorsed BOIN; for borrowing use robust MAP priors; for adaptive platforms use posterior probability of efficacy stopping with simulation-calibrated thresholds.

## Regulatory Status -- The 2024-2026 Bayesian Pivot

**FDA 2010 CDRH Bayesian Devices Guidance** (Feb 5 2010): the only Bayesian-specific FDA guidance until January 2026. Why devices were ahead: CDRH's PMA pathway permits one pivotal trial and accepts borrowing from prior/OUS data more readily than CDER. Example: Edwards SAPIEN (PARTNER B, PMA P100041, Nov 2011) approved using Bayesian propensity-matched comparison to registry of standard-of-care patients.

**FDA January 2026 CDER Bayesian Methodology Draft** (FDA-2025-D-3217; comment period closed March 13 2026): first-ever drug-side Bayesian guidance. Explicit that Bayesian primary inference in pivotals is acceptable provided:
- Prospective specification
- Simulation-based operating characteristics (including frequentist Type-I error under null scenarios — agency still wants calibration)
- Justified priors
- Code/data sufficient for FDA replication

**Project Optimus (FDA OCE, launched 2021; final dose-optimisation guidance Aug 2024):** rewrites Phase I/II oncology by requiring randomised dose comparison before registration. Has made multi-arm randomised dose-finding (BOIN-12, gBOIN-ET) much more important than classic MTD-finding.

**FDA BOIN Fit-for-Purpose qualification (December 2021):** first formal FDA endorsement of a specific dose-finding design under the Drug Development Tools program.

**ICH E20 (Step 2b/3 draft June 2025; NOT final)** treats Bayesian as a legitimate analytic framework but requires demonstration of acceptable frequentist operating characteristics (Type-I, power) over a pre-specified parameter space.

## Algorithmic Taxonomy

| Method | Use case | Software | Strength | Fails when |
|--------|----------|----------|----------|------------|
| BOIN | Phase I MTD | R `BOIN` (Yuan) | **FDA Fit-for-Purpose 2021**; pre-tabulated decisions; no bedside Bayesian software | Statistically less efficient than CRM under correct skeleton |
| mTPI-2 / Keyboard | Phase I MTD | R `escalation`; R `Keyboard` | Default replacement for mTPI; fixes Ockham bias | Tabulated; transparency |
| CRM | Phase I MTD | R `dfcrm`, `trialr` | Most efficient under correct skeleton | Skeleton mis-specification biases MTD |
| EWOC | Phase I MTD | R `ewoc`, `dfcrm` | Explicit overdose-control constraint (P(dose>MTD) <= 0.25) | More conservative than CRM in small trials |
| BOIN-12 / gBOIN-ET | Phase 1b dose-optimisation (Project Optimus) | R `BOIN` extensions | Multi-arm randomised dose comparison | Requires explicit efficacy + toxicity scoring |
| MAP prior | Borrowing from historical control arms | R `RBesT::gMAP` | Industry-standard borrowing | Sample-size of MAP prior must be calibrated (Schmidli 2014) |
| Robust MAP | Borrowing with prior-data conflict protection | R `RBesT::robustify` | Adds vague component (weight 0.1-0.3) to detach if conflict | Mixture weight choice affects borrowing |
| EXNEX | Basket trial across rare-disease strata | R `bhmbasket`; OncoBayes2 | Avoids HM catastrophic borrowing; mixture 0.5/0.5 default (Neuenschwander 2016) | Default weights may over-borrow |
| Dixon-Simon shrinkage | Subgroup analysis | Custom Stan/brms | Honest about no qualitative interaction prior | Prior on tau drives results |
| Berry-Berry 3-level hierarchical | AE multiplicity (AE within PT within SOC) | R `c212`; JMP Clinical | Tames safety multiplicity | Spike-and-slab tuning matters |
| Posterior probability stopping | Adaptive sequential | Custom; FACTS commercial | Bayesian likelihood-principle compatible | Threshold calibration via simulation |
| Predictive probability of success | End-of-Phase-2 go/no-go | Custom Stan | Decision-theoretic; integrates over posterior | Requires Phase 3 design specified |
| Spiegelhalter skeptical/enthusiastic prior | Sensitivity for regulatory pivotals | Custom | Frames regulator-vs-sponsor evidence | Prior elicitation effort |
| Power prior | Pediatric extrapolation borrowing from adults | R `bayesDP`, `psborrow2` | Partial borrowing with discount gamma | gamma choice (Jan 2026 FDA draft: 0.3-0.6) |

**Postdoc reading list:**

- FDA 2010 *Guidance for Industry: Use of Bayesian Statistics in Medical Device Clinical Trials* (Feb 5 2010)
- FDA 2026 Draft *Use of Bayesian Methodology in Clinical Trials* (FDA-2025-D-3217, Jan 2026)
- Berry SM, Carlin BP, Lee JJ, Müller P 2010 *Bayesian Adaptive Methods for Clinical Trials* (CRC)
- Schmidli H, Gsteiger S, Roychoudhury S, O'Hagan A, Spiegelhalter D, Neuenschwander B 2014 *Biometrics* 70:1023 (MAP + robust MAP)
- Weber S, Li Y, Seaman J, Kakizume T, Schmidli H 2021 *J Stat Softw* 100:19 (RBesT)
- Neuenschwander B, Wandel S, Roychoudhury S, Bailey S 2016 *Pharm Stat* 15:123 (EXNEX)
- Liu S, Yuan Y 2015 *J R Stat Soc C* 64:507 (BOIN)
- O'Quigley J, Pepe M, Fisher L 1990 *Biometrics* 46:33 (CRM)
- Babb J, Rogatko A, Zacks S 1998 *Stat Med* 17:1103 (EWOC)
- Ji Y, Liu P, Li Y, Bekele BN 2010 *Clin Trials* 7:653 (mTPI)
- Guo W, Wang SJ, Yang C, Ji Y 2017 *Contemp Clin Trials* 58:23 (mTPI-2 / Keyboard)
- Berry SM, Broglio KR, Groshen S, Berry DA 2013 *Clin Trials* 10:720 (basket trial hierarchical)
- Berry SM, Berry DA 2004 *Biometrics* 60:418 (three-level AE hierarchical)
- Spiegelhalter DJ, Freedman LS, Blackburn PR 1986 *Stat Med* 5:421 (skeptical prior framework)
- Park JW et al 2016 *NEJM* 375:11 (I-SPY 2 veliparib)
- Angus DC et al 2020 *JAMA* (REMAP-CAP COVID rationale)

## Decision Tree by Scenario

| Scenario | Recommended approach | Why |
|----------|---------------------|-----|
| Phase 1 oncology, single-agent MTD | BOIN with target DLT 30%; cohort size 3 | FDA Fit-for-Purpose 2021; tabulated escalation |
| Phase 1 oncology, combination (2 agents) | BLRM with EXNEX in OncoBayes2 | Multi-dimensional dose; industry standard at Novartis/Roche |
| Phase 1b/2 dose-optimisation (Project Optimus) | BOIN-12 or gBOIN-ET; randomised 2-dose comparison | Aug 2024 FDA dose-optimisation guidance |
| Phase 3 with historical control arms available | Robust MAP via RBesT; gMAP() + robustify() | Industry standard borrowing with prior-data conflict protection |
| Basket trial across rare-disease strata | EXNEX (0.5 EX / 0.5 NEX mixture) via OncoBayes2 | Avoids HM catastrophic borrowing |
| Pediatric extrapolation from adult data | Power prior with discount gamma 0.3-0.6 | working convention; the FDA Bayesian Jan 2026 draft does not prescribe a specific gamma range -- check the draft for the current language before quoting |
| Phase 3 trial with single arm + RWE comparator | Propensity-score-integrated power prior via psborrow2 | FDA-supported package for external controls |
| Adaptive trial wanting posterior-probability stopping | Custom Stan model + simulation-calibrated threshold | Bayesian likelihood-principle compatible; no penalty for repeated looks |
| End-of-Phase-2 go/no-go | Predictive probability of success in Phase 3 | Integrates posterior over Phase 3 design |
| Hypothesis-generating safety AE analysis (>100 PTs) | Berry-Berry 3-level hierarchical (AE within PT within SOC) | Tames multiplicity; spike-and-slab on log OR |
| Subgroup analysis post-signal | Bayesian shrinkage (Dixon-Simon, RBesT) | Hemmings-Koch 2019: shrinkage for replication planning, NOT signal generation |
| Regulatory pivotal sensitivity | Spiegelhalter skeptical-prior framework | Frames "evidence for regulators" vs "evidence for sponsor" |

## Phase I Dose-Finding -- BOIN, CRM, mTPI-2

### BOIN (FDA-preferred operational)

```r
library(BOIN)

# Generate escalation table for protocol
boundary_table <- get.boundary(
    target = 0.30,           # target DLT rate
    ncohort = 10,            # 10 cohorts -> max 30 patients with size 3
    cohortsize = 3,
    n.earlystop = 12,        # stop early at lowest dose if 12 patients show futility
    p.saf = 0.6 * 0.30,      # "safe" escalation boundary
    p.tox = 1.4 * 0.30       # "toxic" de-escalation boundary
)
print(boundary_table)
# Pre-printed at investigator desk; no bedside Bayesian software

# Operating characteristics simulation
oc_boin <- get.oc(
    target = 0.30,
    p.true = c(0.05, 0.10, 0.20, 0.30, 0.40, 0.55),  # true DLT per dose
    ncohort = 10,
    cohortsize = 3,
    ntrial = 1000
)
print(oc_boin)
# Reports: MTD selection accuracy, overdose risk, average sample size
```

**BOIN's transparency-over-modelling philosophy:** unlike CRM, BOIN does NOT use information from intermediate dose levels in a model-based way. The Jin-Yuan vs Neuenschwander/Mozgunov debate (Stat Med, Pharm Stat, since ~2018): BLRM/CRM are statistically more efficient under correct model; BOIN is operationally simpler and more transparent.

### CRM with calibrated skeleton

```r
library(dfcrm)

prior_skeleton <- getprior(halfwidth = 0.05, target = 0.30, nu = 3, nlevel = 6)
# Lee-Cheung 2009 indifference-interval calibration

crm_sim <- crmsim(
    PI = c(0.05, 0.10, 0.20, 0.30, 0.40, 0.55),
    prior = prior_skeleton,
    target = 0.30,
    n = 30,
    x0 = 1,                  # starting dose
    nsim = 1000,
    method = 'bayes',
    model = 'logistic'
)
print(crm_sim)
```

**Skeleton mis-specification is the canonical CRM failure mode.** Lee-Cheung 2009 indifference-interval method gives a systematic calibration approach.

### EWOC (overdose control)

```r
# Babb-Rogatko-Zacks 1998: explicit P(dose > MTD) <= alpha (default 0.25)
# Implementation in dfcrm::ewoc; or `ewoc` package
```

## MAP Priors and RBesT

**Schmidli et al 2014 *Biometrics* 70:1023:** Meta-Analytic-Predictive prior. Fit random-effects meta-analysis of historical control arms; derive predictive distribution for new control arm; use as informative prior. Effective sample size from history typically 20-80% of new control arm.

```r
library(RBesT)

# Historical control data (4 prior studies)
historical_data <- data.frame(
    study = c('s1', 's2', 's3', 's4'),
    n = c(40, 35, 50, 45),
    r = c(8, 6, 12, 9)         # responders
)

# Fit MAP via gMAP (Stan-based random-effects meta-analysis)
map_prior <- gMAP(
    cbind(r, n - r) ~ 1 | study,
    data = historical_data,
    family = binomial,
    tau.dist = 'HalfNormal',
    tau.prior = 0.5,           # between-study SD prior
    beta.prior = cbind(0, 2)    # weakly informative on logit response
)
print(map_prior)

# Approximate posterior with mixture for downstream computation
map_mix <- automixfit(map_prior, Nc = 2)
print(map_mix)

# Effective sample size
ess(map_mix)

# Robust MAP: add vague mixture component (weight 0.1-0.3) to guard against prior-data conflict
robust_map <- robustify(map_mix, weight = 0.2, mean = 0.5, n = 1)
print(robust_map)
ess(robust_map)
```

**Robust MAP rationale:** if the new data disagree with historical (prior-data conflict), the mixture down-weights the informative component automatically. Schoenfeld 2017 critique: Schmidli 2014's choice of mixture weight requires sensitivity analysis.

## EXNEX for Basket Trials

**Neuenschwander, Wandel, Roychoudhury, Bailey 2016 *Pharm Stat* 15:123:** Mixture of exchangeable (shared mean+variance) + non-exchangeable (per-basket independent), typically weighted 0.5/0.5. Avoids HM catastrophic borrowing when one basket truly different.

```r
library(OncoBayes2)  # Novartis-developed; canonical EXNEX implementation

# Or simplified via bhmbasket
library(bhmbasket)

# Conceptual: each basket has its own posterior, with shrinkage governed by exchangeability mixture
# Default weights 0.5 EX / 0.5 NEX
# Sensitivity over weights (0.1, 0.3, 0.5, 0.7, 0.9) is essential
```

## Bayesian Platform Trials

### I-SPY 2 (Park-Liu 2016 *NEJM* 375:11)

Neoadjuvant breast cancer; 10 biomarker-defined subtypes × multiple arms; Bayesian RAR; **graduation criterion = posterior predictive probability of success in 300-patient Phase 3 ≥ 0.85.** Berry Consultants designed engine.

```r
# Conceptual implementation requires custom Stan or FACTS (Berry Consultants commercial)

# Pseudocode:
# 1. Fit hierarchical model to platform data: response ~ arm + biomarker_subtype + arm:subtype
# 2. Posterior draws of treatment effect by subtype
# 3. For each draw, simulate Phase 3 trial: n=300, treatment vs control, observed effect
# 4. Compute proportion of draws meeting Phase 3 success criterion
# 5. If proportion >= 0.85, arm graduates
```

### REMAP-CAP (Angus 2020 *JAMA*)

Severe pneumonia, repurposed for COVID-19; Bayesian factorial multi-domain design. Generated corticosteroid signal independently of RECOVERY.

### Drop-the-loser vs promising-the-winner

- **Adaptive arm-dropping (futility):** posterior P(beating control) drops below threshold -> close. Mathematically straightforward.
- **"Promising-the-winner":** selection bias. Bias-adjusted estimators (Robertson 2023; conditional MLE) standard in I-SPY 2 reports.

## Hierarchical Models for Safety Multiplicity (Berry-Berry 2004)

**Berry SM, Berry DA 2004 *Biometrics* 60:418:** three-level hierarchical model for AE multiplicity (AE within MedDRA PT within SOC); spike-and-slab on the log OR. Tames the FDA-feared multiplicity in safety summaries.

```r
library(c212)  # Berry-Berry implementation

# Conceptual: each AE has log OR drawn from spike-and-slab prior
# Spike at 0 (no effect); slab as N(mu_SOC, sigma_SOC)
# SOC-level parameters from N(mu_overall, sigma_overall)
# Borrowing within SOC; shrinkage toward 0 if no evidence

# JMP Clinical also implements this for industry use
```

## Power Priors for Borrowing

```r
library(bayesDP)
library(psborrow2)  # FDA-supported package

# Power prior: combines current data L(theta | D_current) with historical L(theta | D_hist)^gamma
# gamma in [0, 1]; gamma = 0 = no borrowing; gamma = 1 = full pooling

# Typical pediatric extrapolation: gamma = 0.3 to 0.6 per FDA Bayesian Jan 2026 draft
```

## External Control Arms and Real-World Evidence (RWE)

**The 2024-2026 regulatory shift:** FDA has materially expanded acceptance of external/historical/synthetic control arms in rare disease, paediatric, and accelerated-approval settings. Key documents: FDA 2018 RWE Framework (and 2024 enhancements), FDA 2023 Considerations for Use of RWE/RWD for Regulatory Decisions, EMA Reflection Paper on Use of RWE in Regulatory Decision-Making (effective 2024). Bayesian methods are the natural fit because historical data become prior information rather than concurrent control.

### Methodology taxonomy

| Method | Borrowing mechanism | Discount control | When to use |
|--------|---------------------|------------------|-------------|
| Power prior (Ibrahim-Chen 2000) | Likelihood of historical data raised to power gamma | gamma in [0, 1] fixed or modelled | When historical data is single source; gamma ~ Beta in adaptive power prior |
| Robust MAP (Schmidli 2014) | Meta-analytic-predictive prior + vague mixture | Mixture weight (typ 0.1-0.3) | Multiple historical control arms; standard for borrowing |
| Commensurate prior (Hobbs 2011) | Conditional model on agreement parameter | Tau estimated from data | When agreement between historical and current is data-determined |
| Propensity-integrated power prior | Power prior weighted by PS overlap | gamma * (PS-trimmed overlap) | RWE comparator with covariate imbalance |
| Doubly robust ATT via causal inference | IPW + outcome regression | n/a | RWE comparator; identifies marginal ATT |

### psborrow2 — the FDA-supported RWE framework

The `psborrow2` package (Genentech / Bayer / FDA-Janssen collaboration; CRAN 2024+) is the canonical R implementation for propensity-score-integrated Bayesian Dynamic Borrowing. **The skeleton below illustrates the workflow conceptually; verify exact function names and arguments against the current `psborrow2` vignette before use** (the package API has evolved through 2024-2026).

```r
library(psborrow2)

# Define external and internal data
ext_data <- data.frame(usubjid = ..., trt = 0, outcome = ..., covariates = ...)
int_data <- data.frame(usubjid = ..., trt = 0 | 1, outcome = ..., covariates = ...)

# Create borrowing design
borrowing_design <- borrowing_full(
    method_name = "BDB",  # Bayesian Dynamic Borrowing
    ext_flag_col = "ext",
    tau_prior = prior_gamma(0.001, 0.001)  # weakly informative on borrowing
)

# Outcome model (Cox for TTE; logistic for binary)
outcome_model <- outcome_surv_exponential(
    time_var = "time",
    cens_var = "cens",
    baseline_prior = prior_normal(0, 100),
    trt_prior = prior_normal(0, 100)
)

# Run Bayesian analysis with covariate adjustment + borrowing
result <- create_analysis_obj(
    data_matrix = borrow_obj,
    outcome = outcome_model,
    borrowing = borrowing_design,
    covariates = c("age", "ecog", "baseline_severity")
)
mcmc_result <- mcmc_sample(result, n_chains = 4, n_iter = 4000)
```

### Operational rules (FDA 2024-2025 RWE practice)

1. **Pre-specify the RWE source** and document acquisition (registry, EHR, claims, RWD vendor)
2. **Demonstrate comparability** via propensity-score overlap (standardised mean differences <0.25 for key prognostic factors)
3. **Apply discount priors** — full pooling (gamma=1) is regulatory-rejected; typical discount gamma 0.3-0.6
4. **Sensitivity over borrowing strength** — report results at multiple gamma or mixture weights
5. **Tipping-point analysis on prior-data agreement** — at what discount does the conclusion flip?
6. **E-value or bound for unmeasured confounding** (VanderWeele-Ding 2017) — required for FDA submissions; reports the minimum strength of unmeasured confounding that could overturn the result

### When RWE is NOT acceptable

- Trial sponsor and RWE source have meaningful incentive misalignment (e.g., RWE from non-disinterested source)
- RWE captured before standard-of-care evolved (constancy violation, similar to NI biocreep)
- Outcome definitions differ between RWE and current trial (variable harmonisation impossible)
- Censoring patterns in RWE differ structurally from trial (administrative vs disease-driven)
- Highly variable baseline characteristics impossible to balance via propensity weighting

### Recent decisive cases (2024-2026)

- **Zynteglo (FDA 2022, ongoing post-market):** beta-thalassemia gene therapy; single-arm trial vs natural history RWE comparator
- **Skysona (FDA 2022):** cerebral adrenoleukodystrophy; RWE natural-history comparator
- **Multiple ultra-rare disease accelerated approvals 2024-2025:** RWE/external control increasingly accepted in <100-patient trials

## Spiegelhalter Skeptical/Enthusiastic Priors

**Spiegelhalter, Freedman, Blackburn 1986 *Stat Med* 5:421:** the trip-wire / skeptical-prior framework. Pre-specify a skeptical prior centred at the null and an enthusiastic prior centred at the alternative; stopping requires the skeptic to be convinced (posterior under skeptical prior exceeds threshold).

**Frames "evidence for regulators" vs "evidence for sponsor" in Bayesian language**; still cited in modern Bayesian-trial protocols.

```r
# Skeptical prior: N(0, sd_sk) — centred at null
# Enthusiastic prior: N(delta_alt, sd_en) — centred at clinically meaningful effect
# Decision: stop for efficacy if P(theta > 0 | skeptical posterior) > 0.975
#           stop for futility if P(theta < delta_alt | enthusiastic posterior) > 0.80
```

## Per-Method Failure Modes

### CRM with mis-calibrated skeleton

- **Trigger:** Default or arbitrary skeleton without indifference-interval calibration.
- **Mechanism:** Skeleton dictates target dose; mis-calibration biases MTD.
- **Symptom:** MTD selection differs systematically from clinical expectation.
- **Fix:** Calibrate via Lee-Cheung 2009; or switch to BOIN.

### MAP prior with prior-data conflict

- **Trigger:** Historical control rate differs substantially from observed current control.
- **Mechanism:** Informative MAP prior pulls toward historical; current data poorly fit.
- **Symptom:** Posterior dominated by prior; current data evidence under-weighted.
- **Fix:** Robust MAP with mixture weight 0.2-0.3; verify prior-data conflict via posterior predictive checks.

### EXNEX with default 0.5/0.5 weights

- **Trigger:** Default mixture weights without sensitivity.
- **Mechanism:** 50% EX weight allows substantial borrowing even when basket differs.
- **Symptom:** Detected differential basket "softened" by borrowing.
- **Fix:** Sensitivity analysis over weights (0.1, 0.3, 0.5, 0.7, 0.9); report range.

### Posterior probability stopping without simulation-calibrated threshold

- **Trigger:** Stopping rule P(theta > 0 | data) > 0.975 applied without Type-I simulation.
- **Mechanism:** Bayesian rule may not control frequentist Type-I in regulatory sense.
- **Symptom:** FDA review flags lack of Type-I demonstration.
- **Fix:** Simulate under null; calibrate threshold so frequentist Type-I = nominal.

### I-SPY 2 graduation criterion without bias correction

- **Trigger:** Graduated arm's effect estimate reported uncorrected.
- **Mechanism:** Selection on PP > 0.85 inflates estimate.
- **Symptom:** Phase 3 confirmation finds smaller effect than platform suggested.
- **Fix:** Bias-correction via conditional MLE or hierarchical Bayesian; cite Robertson 2023.

### Bayesian shrinkage for signal discovery (Dane vs Hemmings)

- **Trigger:** Hierarchical model fit during signal discovery rather than replication planning.
- **Mechanism:** Shrinkage pre-emptively damps heterogeneity being searched for.
- **Symptom:** Signal detected by causal forest gets shrunken to null in shrinkage analysis.
- **Fix:** Hemmings-Koch 2019 position — shrinkage for replication planning, not signal generation; cite Dane et al 2019 EFSPI white paper + critique.

### Power prior with gamma = 1 (full pooling)

- **Trigger:** Full pooling of historical and current data.
- **Mechanism:** Ignores between-study heterogeneity; biases estimate.
- **Symptom:** Overconfident posterior; cross-validation reveals poor fit.
- **Fix:** Working-convention discount gamma 0.3-0.6 (the FDA Bayesian Jan 2026 draft does not prescribe a specific range); sensitivity over gamma.

### WinBUGS reproducibility

- **Trigger:** Submission contains WinBUGS code without containerised environment.
- **Mechanism:** Older Windows-only software; reproducibility fragile.
- **Symptom:** Reviewer cannot replicate analysis.
- **Fix:** Migrate to Stan (`rstan`/`cmdstanr`); Docker/renv-pinned environment; include seeds + posterior diagnostics (R-hat <1.01, ESS >1000 per chain).

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| FDA BOIN Fit-for-Purpose qualification (Dec 2021) | FDA Drug Development Tools program | First formal FDA dose-finding endorsement |
| Target DLT rate 30% (Phase 1 oncology) | Standard convention | Modal target across oncology Phase 1 |
| MAP prior effective sample size 20-80% of new control | Schmidli 2014 | Borrowing strength typical range |
| Robust MAP mixture weight 0.1-0.3 | Schmidli 2014 | Guards against prior-data conflict |
| EXNEX default 0.5 EX / 0.5 NEX | Neuenschwander 2016 | Standard starting weight; sensitivity required |
| I-SPY 2 graduation PP >= 0.85 | Barker 2009 | Bayesian platform standard |
| Power prior gamma 0.3-0.6 for pediatric extrapolation | working convention; the FDA Bayesian Jan 2026 draft does not prescribe a specific range | Partial borrowing default |
| Stan R-hat <1.01, ESS >1000 per chain | Vehtari 2021 *Bayesian Analysis* | Posterior convergence criteria |
| EWOC overdose constraint P(dose > MTD) <= 0.25 | Babb-Rogatko-Zacks 1998 | Safety floor |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| CRM with arbitrary skeleton | No calibration | Lee-Cheung 2009 indifference-interval; or BOIN |
| MAP without prior-data conflict check | Posterior dominated by prior | Robust MAP; PP-check; sensitivity over mixture weight |
| EXNEX with single weight scheme | No sensitivity | Weights 0.1, 0.3, 0.5, 0.7, 0.9; report range |
| Posterior probability stopping without Type-I sim | Regulatory rejection | Simulate under null; calibrate threshold |
| I-SPY 2 graduated arm reported uncorrected | Selection bias | Conditional MLE; cite Robertson 2023 |
| Bayesian shrinkage for signal discovery | Hemmings-Koch critique | Shrinkage for replication only |
| Power prior gamma = 1 | Full pooling | Discount 0.3-0.6 per FDA 2026 draft |
| WinBUGS without containerisation | Reproducibility | Stan + Docker/renv-pinned |
| BOIN vs CRM comparison without simulation OCs | Apples-to-oranges | Compare OCs over same true DLT rates |
| FDA cited for Bayesian drugs guidance pre-2026 | Confusion | FDA 2010 is DEVICES; FDA 2026 (draft) is drugs |

## Anticipated Reviewer Pushback

| Pushback | Response |
|----------|----------|
| "Type-I error control?" | Simulation under null demonstrates frequentist Type-I = nominal at threshold chosen; documented in SAP appendix |
| "Prior justification?" | MAP from historical control arms via gMAP; robust mixture weight 0.2 for prior-data conflict; sensitivity over prior provided |
| "Why BOIN over CRM?" | BOIN Fit-for-Purpose qualified Dec 2021; pre-tabulated escalation; no bedside Bayesian software; OCs comparable to CRM in simulation |
| "EXNEX weight sensitivity?" | Reported over weights 0.1, 0.3, 0.5, 0.7, 0.9; results stable; primary at 0.5/0.5 per Neuenschwander 2016 |
| "Power prior gamma?" | Discount 0.5 per FDA Bayesian Jan 2026 draft; sensitivity over 0.3-0.7 provided |
| "Posterior probability threshold?" | Calibrated via simulation to frequentist Type-I 0.025 one-sided; cite Berry 2010 |
| "Stan reproducibility?" | Docker container + renv-pinned R + Stan version; seeds provided; R-hat <1.01, ESS >2000 per parameter |
| "Bias correction on platform graduation?" | Conditional MLE applied to estimate Phase 3 effect; cite Robertson 2023 |
| "Why not frequentist instead?" | Bayesian framework permits borrowing (rare disease, pediatric); working convention; the FDA Bayesian Jan 2026 draft does not prescribe a specific gamma range -- check the draft for the current language before quoting primary inference with simulation calibration |

## References

- Babb J, Rogatko A, Zacks S. 1998. Cancer Phase I clinical trials: efficient dose escalation with overdose control. *Stat Med* 17:1103-1120.
- Berry SM, Berry DA. 2004. Accounting for multiplicities in assessing drug safety: a three-level hierarchical mixture model. *Biometrics* 60:418-426.
- Berry SM, Broglio KR, Groshen S, Berry DA. 2013. Bayesian hierarchical modeling of patient subpopulations: efficient designs of Phase II oncology clinical trials. *Clin Trials* 10:720-734.
- Berry SM, Carlin BP, Lee JJ, Müller P. 2010. *Bayesian Adaptive Methods for Clinical Trials*. CRC.
- Dane A, Spencer A, Rosenkranz G, Lipkovich I, Parke T. 2019. Subgroup analysis and interpretation for phase 3 confirmatory trials: EFSPI/PSI white paper. *Pharm Stat* 18:126-139.
- FDA. 2010. Guidance for Industry: Use of Bayesian Statistics in Medical Device Clinical Trials.
- FDA. 2021. BOIN Drug Development Tool Fit-for-Purpose Qualification.
- FDA. 2026. Use of Bayesian Methodology in Clinical Trials. Draft Guidance (FDA-2025-D-3217).
- Guo W, Wang SJ, Yang C, Ji Y. 2017. A Bayesian interval dose-finding design addressing Ockham's razor: mTPI-2. *Contemp Clin Trials* 58:23-33.
- Hemmings R, Koch A. 2019. Commentary on Dane et al. *Pharm Stat* 18:140-141.
- Ji Y, Liu P, Li Y, Bekele BN. 2010. A modified toxicity probability interval method for dose-finding trials. *Clin Trials* 7:653-663.
- Liu S, Yuan Y. 2015. Bayesian optimal interval designs for phase I clinical trials. *JRSS-C* 64:507-523.
- Neuenschwander B, Wandel S, Roychoudhury S, Bailey S. 2016. Robust exchangeability designs for early phase clinical trials with multiple strata. *Pharm Stat* 15:123-134.
- O'Quigley J, Pepe M, Fisher L. 1990. Continual reassessment method: a practical design for phase 1 clinical trials in cancer. *Biometrics* 46:33-48.
- Park JW et al. 2016. Adaptive randomization of veliparib-carboplatin treatment in breast cancer. *NEJM* 375:11-22.
- Robertson DS, Lee KM, López-Kolkovska BC, Villar SS. 2023. Response-adaptive randomization in clinical trials: from myths to practical considerations. *Stat Sci* 38:185-208.
- Schmidli H, Gsteiger S, Roychoudhury S, O'Hagan A, Spiegelhalter D, Neuenschwander B. 2014. Robust meta-analytic-predictive priors in clinical trials with historical control information. *Biometrics* 70:1023-1032.
- Spiegelhalter DJ, Freedman LS, Blackburn PR. 1986. Monitoring clinical trials: conditional or predictive power? *Stat Med* 5:421-433.
- Vehtari A et al. 2021. Rank-normalization, folding, and localization: an improved R-hat for assessing convergence. *Bayesian Analysis*.
- Weber S, Li Y, Seaman J, Kakizume T, Schmidli H. 2021. Applying meta-analytic-predictive priors with the R Bayesian evidence synthesis tools. *J Stat Softw* 100:19.

## Related Skills

- clinical-biostatistics/adaptive-designs - Group-sequential, SSR, platform trials
- clinical-biostatistics/subgroup-analysis - Bayesian shrinkage for HTE (Dixon-Simon, Berry)
- clinical-biostatistics/power-and-sample-size - Bayesian SS via predictive probability of success
- clinical-biostatistics/multiplicity-graphical - Berry-Berry AE hierarchical
- clinical-biostatistics/trial-reporting - Bayesian inference reporting per CONSORT 2025
- clinical-biostatistics/missing-data-sensitivity - Bayesian rbmi imputation
- machine-learning/biomarker-discovery - Bayesian HTE for biomarker subgroups
- experimental-design/sample-size - General methods
