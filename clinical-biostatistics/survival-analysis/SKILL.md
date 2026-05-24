---
name: bio-clinical-biostatistics-survival-analysis
description: Performs time-to-event analysis for clinical trials including Cox proportional hazards regression with PH diagnostics, restricted mean survival time (RMST) under non-PH, competing risks via Fine-Gray vs cause-specific Cox, weighted log-rank and MaxCombo for non-proportional hazards, recurrent events (Andersen-Gill, PWP, WLW), and interval-censored data. Use when analyzing time-to-event endpoints (OS, PFS, DOR, TTR, TTNT) in oncology or other clinical trials.
tool_type: mixed
primary_tool: lifelines
---

## Version Compatibility

Reference examples tested with: lifelines 0.27+, scikit-survival 0.21+, statsmodels 0.14+, pandas 2.1+, numpy 1.26+. R packages cited (still the SOTA for survival): survival 3.8+, survRM2, cmprsk, riskRegression, mstate, flexsurv, icenReg, rpsftm.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Time-to-Event Analysis for Clinical Trials

**"Analyze time-to-event endpoint"** -> Estimate a hazard, survival probability, cumulative incidence, or restricted mean time using a method calibrated to (a) whether proportional hazards holds, (b) whether competing events exist, (c) whether censoring is informative, and (d) which estimand the trial targets under ICH E9(R1).

## The Single Most Important Modern Insight -- PH Almost Never Holds

In modern oncology with checkpoint inhibitors, targeted therapies, crossover, and depleted high-risk subjects over follow-up, **proportional hazards (PH) violations are the rule, not the exception**. The Cross-Pharma NPH Working Group (Lin et al 2020 *Stat Biopharm Res*; Magirr-Burman 2021 *Stat Biopharm Res* 15(2):293) documented systematic PH violations across phase III oncology trials, particularly delayed-effect patterns from checkpoint inhibitors.

The Cox HR is a *time-averaged log-hazard ratio* under PH violation (Xu-O'Quigley 2000), which may or may not be the estimand of interest. RMST (Royston-Parmar 2013) provides a clinically interpretable, hazard-free alternative.

## Algorithmic Taxonomy

| Method | Estimand | Inference | Strength | Fails when |
|--------|----------|-----------|----------|------------|
| Log-rank (unstratified) | Test of S_A(t) = S_B(t) all t | Permutation / asymptotic chi-square | Standard; preserves Type-I under PH | Underpowered under non-PH; treats all events equally |
| Stratified log-rank | Same null within strata, pooled | Asymptotic | Preserves stratification factor from randomisation | Stratification factor must be pre-specified |
| Weighted log-rank G(rho, gamma) | Direction-specific test under non-PH | Asymptotic | High power for delayed/early/middle effects | Weight choice must match true effect time profile; chasing weight = p-hacking |
| Cox PH | Conditional log-HR | Wald, LR, score | Standard; semi-parametric; covariate adjustment | PH violation makes HR a misleading summary; check via cox.zph |
| Stratified Cox | Same; baseline hazards differ by stratum | Wald | Handles non-PH by stratification | Loses inference on stratification variable; cannot interact treatment with strata |
| Time-varying Cox (`tt()`) | Time-dependent log-HR | Wald | Quantifies non-PH explicitly | Interpretability — no single "the HR"; choose g(t) carefully |
| Flexible parametric (Royston-Parmar) | Time-varying log-HR via splines | Wald | Smooth S(t), HR(t); supports extrapolation | Spline choice affects results; software in R `stpm2/stpm3` |
| RMST | Difference in mean survival truncated at tau | Wald with delta or pseudo-obs regression | Hazard-free; clinically interpretable in time units | tau choice; min follow-up across arms constrains tau |
| MaxCombo | Maximum over weighted log-rank family | Asymptotic multivariate normal | Robust to range of NPH patterns | Can reject in opposite directions on same data (Magirr 2022 critique) |
| Fine-Gray subdistribution HR | Conditional subdistribution HR | Wald | Direct CIF modeling | Andersen-Keiding 2012 critique: violates causal hazard semantics |
| Cause-specific Cox | Conditional cause-specific HR | Wald | Causally interpretable | Predicts hazards, not CIFs; need both for CIF prediction |
| Multi-state Cox (mstate) | Transition-specific HRs | Wald | Subsumes competing risks; handles relapse/remission | More complex; more parameters to estimate |
| Andersen-Gill (recurrent) | Rate ratio | Robust (cluster) Wald | Most efficient under exchangeability | Assumes exchangeable events |
| PWP (recurrent) | Conditional event-order HR | Stratified Wald | Handles event-order qualitative heterogeneity | More strata = more parameters; smaller per-stratum n |
| Interval-censored Cox (NPMLE) | Cumulative hazard | Likelihood ratio | Correct for periodic-assessment data | Slower; software in R `icenReg` |

**Postdoc reading list:**

- Royston-Parmar 2013 *BMC Med Res Methodol* 13:152 (RMST as primary)
- Uno et al 2014 *JCO* 32:2380 (RMST in oncology)
- Andersen-Keiding 2012 *Stat Med* 31:1074 (Fine-Gray semantic critique)
- Putter-Schumacher-van Houwelingen 2020 *Biom J* 62:790 (Fine-Gray revisited; reduction factor)
- Putter-Fiocco-Geskus 2007 *Stat Med* 26:2389 (competing risks tutorial)
- Magirr-Burman 2021 *Stat Biopharm Res* 15(2):293 (MaxCombo critique)
- Grambsch-Therneau 1994 *Biometrika* 81:515 (scaled Schoenfeld residuals)
- Buyse-Molenberghs 1998 *Biometrics* 54:1014 (PFS-OS surrogacy framework)
- Lewis et al 2023 (ICH E9(R1) censoring rules and estimand)
- Sun 2006 *The Statistical Analysis of Interval-Censored Failure Time Data* (Springer)

## Decision Tree by Scenario

| Scenario | Recommended approach | Why |
|----------|---------------------|-----|
| OS, drug expected to extend survival uniformly, PH plausible | Stratified log-rank + Cox HR with cox.zph diagnostic | Standard regulatory approach; pre-specify stratification factors from randomisation |
| PFS in immuno-oncology with expected delayed separation | MaxCombo with pre-specified G(0,0), G(0,1), G(1,0), G(1,1); RMST as sensitivity | Robust to NPH; explicit direction check (Magirr-Burman 2021) |
| OS with crossover to active arm | Treatment policy estimand (ITT) primary; RPSFT/IPCW as sensitivity for hypothetical | Lewis 2023; FDA/EMA accept both, ITT is primary |
| Time-to-event with assessment-schedule artifact (periodic scans) | Interval-censored Cox (R `icenReg::ic_par`) | Standard right-censoring at midpoint is biased |
| DOR (duration of response) | KM among responders; censor at next-therapy/death/dropout per pre-specified estimand | Weber 2023 *Pharm Stat*; responder-conditioned = doubly post-randomisation |
| Competing risk: drug efficacy on event A in context of high event-B mortality | Cause-specific Cox for A AND for B; report both; CIF via Aalen-Johansen | Andersen-Keiding 2012; Fine-Gray SHR is not causal |
| Prediction of CIF for clinical decision | Fine-Gray for CIF estimate; cause-specific Cox for etiology | Putter-Fiocco-Geskus 2007 split |
| Multi-state model (alive -> relapse -> death) | mstate framework with transition-specific Cox models | Subsumes competing risks |
| Recurrent events (exacerbations, hospitalisations) | Andersen-Gill with robust SE; cite event-order assumption | Most efficient; falls back to PWP if order matters |
| Single-arm OS extrapolation for HTA | Flexible parametric (Royston-Parmar) + external information | Smooth tail; supports extrapolation beyond trial follow-up |
| Non-PH with crossing hazards | RMST with pre-specified tau; OR multi-state model | HR loses meaning under crossing |

## Cox PH Diagnostics -- The Therneau-Grambsch Test (and Its Pitfalls)

**Goal:** Detect violations of the proportional hazards assumption that would invalidate the Cox HR as a meaningful summary statistic.

**Approach:** Compute scaled Schoenfeld residuals (Grambsch-Therneau 1994); regress against a time-transform g(t) under H0 of zero slope; supplement the asymptotic p-value with a graphical residual plot since the test is sensitive to g(t) choice and sample-size dependent.

```python
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test

cph = CoxPHFitter()
cph.fit(df, duration_col='time', event_col='event', formula='treatment + age + baseline_score')
cph.print_summary()

# Schoenfeld residuals PH test (lifelines)
results = proportional_hazard_test(cph, df, time_transform='rank')
print(results.summary)

# In R: cox.zph(coxph_fit) returns same with KM transform as default
```

**The g(t) choice trap (Park-Hendry 2015 *Political Analysis*):** the test power depends critically on g(t). KM transform, identity, log, and rank give materially different p-values. With n > 5000 even trivial deviations reject; with n < 100 the test misses meaningful violations.

**Critical interpretation rule: a global p > 0.05 does NOT mean PH holds — it means the null cannot be rejected.** The graphical diagnostic is more informative — a flat smoothed line is the target. Use cox.zph as a *failure detector*, not a *PH validator*.

```python
# Plot scaled Schoenfeld residuals (lifelines)
cph.check_assumptions(df, p_value_threshold=0.05, show_plots=True)
```

### Fixes when PH violated

1. **Stratify** on the violating covariate (loses inference on it)
2. **Time-dependent coefficients:** `coxph(Surv(t, d) ~ x + tt(x), tt=function(x,t,...) x*log(t))` in R; lifelines: `formula='treatment * time'`
3. **Royston-Parmar flexible parametric** (R `flexsurv::flexsurvspline` or `rstpm2::stpm2`)
4. **RMST** as primary; HR as secondary

## Restricted Mean Survival Time (RMST) -- The Modern Alternative

**Royston-Parmar 2013 + Uno 2014 case:** RMST(tau) = E[min(T, tau)] = integral from 0 to tau of S(t) dt = area under KM curve up to tau. The difference RMST_A(tau) - RMST_B(tau) is a **time gained in months** -- clinically interpretable, requires no PH assumption, always estimable up to the minimum of the largest follow-up across arms.

```python
# Python: rmst via lifelines or manual
from lifelines.utils import restricted_mean_survival_time
from lifelines import KaplanMeierFitter

kmf_A = KaplanMeierFitter().fit(df[df['arm']=='A']['time'], event_observed=df[df['arm']=='A']['event'])
kmf_B = KaplanMeierFitter().fit(df[df['arm']=='B']['time'], event_observed=df[df['arm']=='B']['event'])
tau = 36  # months; pre-specified
rmst_A = restricted_mean_survival_time(kmf_A, t=tau)
rmst_B = restricted_mean_survival_time(kmf_B, t=tau)
print(f'RMST diff: {rmst_A - rmst_B:.2f} months')

# R: survRM2::rmst2(time, status, arm, tau=36) is the standard
```

### The tau (truncation time) choice

- **Statistical constraint:** tau <= min(largest follow-up time in each arm) to avoid extrapolation (Tian et al 2020 *Biostatistics* gives rigorous treatment)
- **Clinical constraint:** tau should reflect a clinically meaningful horizon (5-year OS in adjuvant; 24-month PFS in metastatic)
- **Data-dependent tau inflates Type-I error** — MUST be pre-specified in SAP. Post-hoc tau tuning to chase significance is p-hacking.

### Pseudo-observation regression (Andersen-Klein 2007)

For each subject i, compute jackknife pseudo-value θ_i(tau) = n·RMST(tau) - (n-1)·RMST_{-i}(tau). Regress pseudo-values on covariates via GEE — enables RMST regression WITH covariate adjustment, including time-varying covariates, no PH assumption.

R implementation: `pseudo::pseudomean()` + GEE via `geepack::geeglm()`.

**Lambert critique (postdoc reading):** Paul Lambert (`stpm3` author) is on record that RMST is "not unambiguously simpler than HR for clinicians." HR has 50-year pedagogy lead. RMST is cumulative (not instantaneous) and cannot detect crossing hazards beyond chosen tau. The argument for RMST is hazard-free interpretability under non-PH; the argument against is communication and tau sensitivity.

## Competing Risks -- The Andersen-Keiding Framework

### Fine-Gray subdistribution hazard

```python
# Python: scikit-survival or lifelines (limited support); R is the SOTA
# R: library(cmprsk); crr(time, fstatus, cov, failcode=1, cencode=0)
# R: library(riskRegression); FGR(formula, data=df, cause=1)
```

**Fine-Gray 1999** introduced subdistribution hazard lambda^FG(t) = -d/dt log[1 - F_1(t)] where F_1 is CIF for cause 1. The proportional subdistribution hazards model places covariates on this hazard. Subjects who experience competing events remain in the risk set with IPCW weighting — mathematically convenient but **conceptually weird**.

**The Andersen-Keiding 2012 critique** (*Stat Med* 31:1074): Fine-Gray subdistribution hazard violates the three principles for valid hazard functionals:

1. Hazard must be a real instantaneous risk among those truly at risk -> FG keeps dead competing-event subjects "at risk"
2. Covariate effects must be interpretable as causal -> FG coefficients confound the competing-event hazard
3. Hazard must support landmarking (conditioning on survival to s) -> FG does not

**The competing-risk confounding trap** (Putter-Schumacher-van Houwelingen 2020 *Biom J* 62:790; the reduction factor decomposition makes it explicit): any covariate that increases the cause-specific hazard of competing event A will *decrease* the subdistribution hazard for event B simply because A removes subjects from the population at risk of B. A Fine-Gray "protective" effect on B may be an iatrogenic killer via A.

### Practical rule (Putter-Fiocco-Geskus 2007)

- **For prediction of CIF (cumulative incidence)** -> Fine-Gray is fine; report CIF curves and SHR
- **For etiology / causal effect on the event** -> use cause-specific Cox (treat competing events as censoring); report cause-specific HRs for BOTH event of interest AND competing events
- **For multi-state semantics** (alive -> relapse -> death) -> use multi-state framework with transition-specific Cox models

### The CIF is always estimable

Aalen-Johansen estimator (multi-state generalisation of KM) gives the CIF non-parametrically. **The KM estimator is biased upward in the presence of competing risks** — it treats competing events as non-informative censoring and overestimates 1 - CIF.

```r
# R: library(mstate); msfit(coxph_fit, newdata, trans); probtrans(msfit_obj)
# R: library(survival); survfit(Surv(time, event_factor) ~ 1) with multi-state Surv
```

## Log-Rank Variants and MaxCombo

```python
from lifelines.statistics import logrank_test, multivariate_logrank_test

# Stratified log-rank (lifelines uses `weightings=` -- note plural)
results = multivariate_logrank_test(
    df['time'], df['arm'], df['event'], weightings='peto'
)
# weightings options (lifelines 0.27+): None, 'wilcoxon', 'tarone-ware',
# 'peto', 'fleming-harrington' (also accepts kwargs for Fleming-Harrington p, q)
```

| Test | Citation | When to use |
|------|----------|-------------|
| Standard log-rank (G(0,0)) | Mantel 1966 | PH holds |
| Wilcoxon | Wilcoxon 1945 | Down-weights late events; early-effect detection |
| Tarone-Ware | Tarone-Ware 1977 | Compromise |
| Peto-Peto | Peto-Peto 1972, Prentice 1978 | Robust to ties and censoring distribution differences between arms |
| Fleming-Harrington G(rho, gamma) | Fleming-Harrington 1981 | Direction-specific: G(0,1) for late-emphasis (delayed effects), G(1,0) for early |
| MaxCombo | Karrison 2016 | Take maximum over family; control multiplicity via joint MVN |

### The MaxCombo controversy

**Magirr-Burman 2021 *Stat Biopharm Res* 15(2):293**: MaxCombo can reject the null in **opposite directions** on the same dataset. Formally, it rejects the strong null H_0: S_A(t) = S_B(t) for all t, but the rejection direction is determined by the dominant weight, which can flip across portions of the curve. **KEYNOTE-042 demonstration**: MaxCombo simultaneously favouring pembrolizumab AND chemo depending on weight choice.

**Cross-Pharma NPH Working Group recommendation:** MaxCombo with directionality constraints — require positive z-statistic at the late-emphasis weight before declaring superiority; report the dominant weight and its direction.

## Interval Censoring -- When Standard Right-Censoring Is Wrong

**PFS is REALLY interval-censored** — events known only between consecutive RECIST scans. The convention is to treat it as right-censored at midpoint or first-PD-scan date. This is a long-standing methodological compromise defensible only when scan intervals are short and balanced across arms.

```r
# R: library(icenReg)
# fit_ic <- ic_par(cbind(left, right) ~ treatment + age, data=df, dist='weibull', model='ph')
# ic_sp() for semi-parametric; ic_npar() for NPMLE
```

**Sun 2006** *The Statistical Analysis of Interval-Censored Failure Time Data* (Springer) is the canonical reference. NPMLE via Turnbull's algorithm; Cox-like regression via Finkelstein 1986 / Pan 1999 EM algorithms.

**When to switch from right-censored midpoint to interval-censored:** scan intervals > 4 weeks; scan timing differs between arms (open-label trials with potential differential ascertainment); regulatory submission where 2-3% effect-size shift matters.

## Recurrent Events -- AG vs PWP vs WLW

| Model | Risk set | Baseline hazard | Best when |
|-------|----------|-----------------|-----------|
| Andersen-Gill 1982 | Total time; subject at risk continuously between events | Common (counting process) | Events exchangeable; rate model; no event-order effect |
| PWP gap-time (Prentice-Williams-Peterson 1981) | Stratified by event number; subject enters stratum k after event k-1 | Stratum-specific | Event-order matters; later events qualitatively different |
| WLW (Wei-Lin-Weissfeld 1989) | Marginal -- separate Cox per event order | Per-event-order | Multiple types of events; uses sandwich variance |

```python
# Python: lifelines does not have native AG; use coxph with (start, stop, event) format
# R: coxph(Surv(start, stop, event) ~ trt + cluster(id), data=long_df) for AG
# R: + strata(enum) for PWP
```

**Box-Steffensmeier critique:** WLW is often misused — fitting separate Cox to "time to 2nd event" treats it as a first-event problem rather than conditional on prior history, inflating effect estimates. PWP-gap-time is the cleanest conditional model. AG most efficient when exchangeability holds.

**Modern advice (Rogers et al 2016; Cook-Lawless 2007 book):** AG with robust variance as default; PWP if events qualitatively heterogeneous; avoid WLW unless events are truly distinct types.

## Oncology PFS/OS Estimands -- The 2024 Reality

Per ICH E9(R1), the same PFS dataset yields different HRs depending on censoring rules. Each rule corresponds to a different intercurrent-event strategy:

| Censoring rule | ICE strategy | Estimand |
|----------------|-------------|----------|
| Censor at last assessment before missing visits | Hypothetical (had visit not been missed) | What would PFS be if visits never missed? |
| Censor at start of new anticancer therapy | Hypothetical (no subsequent therapy) | What would PFS be without rescue? |
| Count subsequent therapy as event | Composite (subsequent therapy = treatment failure) | Treatment-failure-free survival |
| No censoring (event = last assessment + progression) | Treatment policy | Real-world PFS with policy of allowing rescue |

**The 2024 European Journal of Cancer demonstration (PMID 38547775):** two sets of censoring rules — FDA-favoured vs trialist-favoured — applied to the *same* PFS data shifted median PFS from 32 to 43 months in the experimental arm with no change in control. This is the estimand changing, not analytic artefact.

**Fleming 2025 argument:** PFS should be defined as treatment-policy with progression + death + subsequent therapy as a composite ("treatment failure-free survival") to preserve ITT. Controversial because it loses interpretation as "tumor growth control."

### Informative censoring -- detection and handling

**The mechanism:** standard right-censored Cox / KM assumes censoring is non-informative (the censoring process is unrelated to the underlying failure time after conditioning on observed covariates). **When this fails, KM is biased upward** in the arm with informative censoring -- patients who drop out due to lack of efficacy or toxicity are precisely the ones who would have failed early.

**TROPiCS-02 2023 (Tolaney *JCO*):** "evaporative cooling" of progression events. Patients on the toxic arm discontinue and are censored BEFORE the next protocol-mandated scan that would have captured progression; KM biases PFS upward in the toxic arm. Templeton 2020 *Nat Rev Clin Oncol* and Campigotto-Weller 2014 *JCO* are foundational citations.

**Detection workflow:**

1. **Tabulate censoring reasons by arm** from CDISC DS (Disposition) or ADaM ADTTE CNSR integer values:
   - CNSR=1 (lost to follow-up): symmetric across arms = OK
   - CNSR=2 (withdrew consent): if asymmetric, investigate why
   - CNSR=3 (admin EoS): symmetric by definition
   - CNSR=4 (subsequent therapy initiated): differential is the canonical informative pattern
   - CNSR=5+ (toxicity discontinuation): differential is highly informative

2. **Compare KM curves for censoring distribution** by arm. If "time-to-censoring" KM differs by arm in same direction as outcome, suspect informative censoring.

3. **Cox-Snell or Schoenfeld residuals on censoring-as-event model**: if treatment is significantly associated with censoring hazard, censoring is informative.

**Handling strategies (per ICH E9(R1) and Lewis 2023):**

| Strategy | When to use | Implementation |
|----------|-------------|----------------|
| Composite endpoint (treatment-policy) | ICE has clinical meaning (subsequent therapy = treatment failure) | Re-define event as "first of (progression OR subsequent therapy OR death)"; eliminates need for censoring assumption |
| IPCW with stabilised weights (Robins 1992) | Hypothetical estimand under "no informative dropout" | R `ipw::ipwtm` or custom: weight = inverse probability of remaining uncensored given baseline + time-varying covariates |
| Sensitivity: worst-case / best-case imputation | Bounding the true effect under informative censoring | Censored patients = events at censoring date (worst case); or remain at risk (best case); report range |
| Multistate model | Multiple competing causes of censoring | `mstate` for transition-specific Cox; handles death + dropout + subsequent therapy as separate states |
| RPSFT / structural nested failure time | Crossover-induced informative censoring | R `rpsftm` package; FDA / EMA acceptable for OS hypothetical estimand under crossover |

**IPCW workflow (stabilised weights per Robins-Finkelstein 2000):**

```r
# Stabilised IPCW: w(t) = S_num(t) / S_denom(t)
# where S_denom predicts remaining-uncensored given baseline + time-varying covariates,
# and S_num predicts remaining-uncensored given treatment alone (numerator for stabilisation)

# Step 1: fit censoring (NOT event) hazard models
# Note: outcome variable for these models is "1 if censored, 0 if event or still at risk"
denom_cens <- coxph(Surv(time, censored) ~ treatment + age + baseline_severity, data=df)
num_cens   <- coxph(Surv(time, censored) ~ treatment, data=df)

# Step 2: compute survival probabilities (NOT hazards) for "remaining uncensored"
# basehaz + linear predictor -> S(t) per subject; the `survfit` + `summary` API gives S_i(t_i)
S_denom <- summary(survfit(denom_cens, newdata=df), times=df$time)$surv
S_num   <- summary(survfit(num_cens,   newdata=df), times=df$time)$surv

# Step 3: stabilised IPCW weight per subject at their observed time
df$ipcw <- S_num / S_denom

# Step 4: weighted Cox for the event of interest with robust SE
fit_ipcw <- coxph(Surv(time, event) ~ treatment + age + baseline_severity,
                  data=df, weights=df$ipcw, robust=TRUE)
```

**Note:** in practice use the `ipw::ipwtm` (or `ipwExt`) wrapper which handles the time-varying weights and stabilisation correctly; the pseudocode above shows the underlying machinery. Verify weight distribution: median weight near 1, no extreme values (>10 indicates near-violation of positivity).

**Operational rule:** for any TTE primary analysis where DS shows differential censoring reasons by arm, the regulatory expectation is (1) tabulate reasons in CSR, (2) report sensitivity under IPCW or composite, (3) discuss whether primary HR/RMST changes under sensitivity.

## ADaM ADTTE -- The CNSR Convention Trap

**CDISC ADTTE convention: CNSR = 0 for events, positive integers for censoring reasons** (1 = lost to follow-up, 2 = withdrew consent, 3 = admin EoS, 4 = subsequent therapy, etc.). **Opposite to most stat packages** which use 1 = event.

```python
# Convert ADTTE for Python / R
import pandas as pd

adtte = pd.read_csv('ADTTE.csv')
adtte['event'] = (adtte['CNSR'] == 0).astype(int)  # 1 = event for survival packages
adtte['time'] = adtte['AVAL']  # AVAL is in days/months per AVALU
```

This is a perpetual bug source. See clinical-biostatistics/cdisc-data-handling for the full ADTTE specification.

## Per-Method Failure Modes

### Cox PH violation undetected

- **Trigger:** Significant treatment HR reported without cox.zph diagnostic
- **Mechanism:** Cox HR is a time-averaged log-HR under PH violation; the "the HR" interpretation breaks
- **Symptom:** Hazard plots show crossing; cox.zph rejects PH; KM curves diverge then converge
- **Fix:** Report cox.zph result; switch to RMST or time-varying Cox; cite Grambsch-Therneau 1994

### Fine-Gray reported as causal effect

- **Trigger:** Fine-Gray SHR interpreted as "treatment effect on event of interest"
- **Mechanism:** SHR confounds with competing-event hazard (Andersen-Keiding 2012)
- **Symptom:** SHR interpretation contradicts cause-specific Cox results; reviewer confusion
- **Fix:** Use cause-specific Cox for etiology; Fine-Gray for CIF prediction only; report both per Putter-Fiocco-Geskus 2007

### KM curve biased upward under competing risks

- **Trigger:** KM applied to event of interest treating competing events as non-informative censoring
- **Mechanism:** KM estimates 1 - cause-specific hazard cumulative; competing events deplete denominator
- **Symptom:** KM survival exceeds Aalen-Johansen 1 - CIF estimate
- **Fix:** Use Aalen-Johansen estimator for CIF; never KM in competing-risk setting

### MaxCombo direction flipping

- **Trigger:** MaxCombo significant; direction not pre-specified
- **Mechanism:** MaxCombo's family can include weights that favour opposite directions
- **Symptom:** Different reports show opposite "winner" depending on weight choice
- **Fix:** Pre-specify directional constraints; require positive z at late-emphasis weight before claiming superiority; cite Magirr-Burman 2021

### Tau chosen post-hoc to favour significance

- **Trigger:** RMST tau adjusted after seeing data
- **Mechanism:** Selection bias; tau-tuning is a flavour of p-hacking
- **Symptom:** Sponsor's RMST result differs from independently re-analysed with pre-specified tau
- **Fix:** Pre-specify tau in SAP; cite Tian 2020 *Biostatistics*

### Right-censored midpoint for interval-censored PFS

- **Trigger:** Scan intervals long or differential between arms; standard right-censored Cox applied
- **Mechanism:** Midpoint imputation is biased; SE underestimates
- **Symptom:** Replication with interval-censored Cox gives different point estimate and wider CI
- **Fix:** R `icenReg::ic_par` or `ic_sp`; Sun 2006

### CNSR convention confusion

- **Trigger:** ADTTE CNSR=1 passed to R `survival::Surv(time, event)` expecting event=1
- **Mechanism:** Censoring/event role reversed
- **Symptom:** "Event count" matches censoring count from CSR; nonsensical HR
- **Fix:** Always convert: `event = (CNSR == 0).astype(int)`; cite ADaM IG v1.3

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Cox HR significant, RMST difference non-significant | HR is time-averaged log-HR under PH violation; RMST captures cumulative effect | If PH violated, RMST is the more interpretable summary; HR may be artifact of late-event preponderance |
| MaxCombo rejects null but direction depends on weight | MaxCombo's rejection direction is determined by dominant weight; can flip across weight choices (Magirr-Burman 2021) | Pre-specify directional constraint (positive z at late-emphasis weight) before declaring superiority; do NOT post-hoc select winning weight |
| Fine-Gray SHR vs cause-specific Cox HR conflict on direction | Fine-Gray confounds with competing-event hazard; cause-specific isolates causal effect on event of interest (Andersen-Keiding 2012) | For causal claims use cause-specific Cox; Fine-Gray for CIF prediction only; report both per Putter-Fiocco-Geskus 2007 |
| Stratified log-rank p < unstratified p | Stratification removes between-stratum variance; correct standard error smaller | Stratified analysis matches randomisation; unstratified inflates Type-I per Kahan-Morris 2012 |
| KM survival curve vs Aalen-Johansen 1-CIF disagree in competing-risks setting | KM treats competing events as non-informative censoring (biased upward) | Use Aalen-Johansen for CIF; never KM in competing-risk setting |
| Schoenfeld-formula SS insufficient at trial end (events under-collected) | Non-PH (immuno-oncology delayed effect) violates Schoenfeld assumption; under-estimates events 20-50% | Re-power using Lakatos 1988 or simulation under expected HR(t); cite Lin 2020 NPH Working Group |
| Right-censored midpoint Cox HR vs interval-censored Cox HR differ | Midpoint imputation biased when scan intervals long or asymmetric | Switch to interval-censored Cox via icenReg when scan intervals > 4 weeks or differ by arm |
| ADTTE CNSR=1 produces nonsensical results | CDISC convention reversal: CNSR=0 means event in ADaM | Always convert: `event = (CNSR == 0).astype(int)` before stat-package call |

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| cox.zph p > 0.05 is NOT proof of PH | Grambsch-Therneau 1994; Park-Hendry 2015 | Failure detector, not validator |
| RMST tau <= min(largest follow-up per arm) | Tian 2020 *Biostatistics* | Avoid extrapolation; tau pre-specified in SAP |
| MaxCombo with directional constraint | Magirr-Burman 2021 *Stat Biopharm Res* 15(2):293 | Prevents opposite-direction rejections |
| Cause-specific Cox for etiology; FG for CIF prediction | Putter-Fiocco-Geskus 2007 | Andersen-Keiding 2012 critique |
| Scan intervals > 4 weeks -> interval-censored analysis | Sun 2006 | Midpoint right-censoring is biased |
| 10 events per covariate for Cox | Peduzzi 1995 *J Clin Epidemiol* | Below this, bias and overfitting |
| Schoenfeld 1981 events formula assumes PH | Schoenfeld 1981 | Under non-PH, under-estimates required events by 20-50%; cite Lakatos 1988 for non-PH SS |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Cox HR reported with cox.zph p < 0.001 ignored | Diagnostic skipped | Report RMST or time-varying Cox; cite Therneau-Grambsch |
| KM curve labeled "survival from event of interest" with competing risks | Bias upward | Aalen-Johansen CIF; never KM with competing risks |
| Fine-Gray HR reported as "treatment effect" | Misinterpretation of SHR | Cite Andersen-Keiding; report cause-specific too |
| MaxCombo with no direction restriction | Direction-flipping risk | Pre-specify constraints (Magirr-Burman 2021) |
| ADTTE CNSR=1 used as "1=event" | Convention reversal | `event = (CNSR == 0).astype(int)` |
| PFS midpoint right-censored without interval-censored sensitivity | Scan-schedule bias | Interval-censored analysis when scan intervals long |
| Stratified randomisation factor not in Cox | Type-I inflation | Include strata via `strata()` or as covariate |
| Schoenfeld SS calculation in immuno-oncology | PH assumption violated | Simulate under expected hazard pattern (Lakatos 1988); use MaxCombo SS |
| tau set to longest follow-up post-hoc | RMST p-hacking | Pre-specify tau in SAP; cite Tian 2020 |
| Recurrent-event WLW with naive "time to 2nd event" | Inflated effect | AG with robust SE, or PWP gap-time |

## Anticipated Reviewer Pushback

| Pushback | Response |
|----------|----------|
| "PH assumption check?" | cox.zph result reported with caveat that p>0.05 doesn't prove PH; graphical diagnostic in supplement |
| "Why RMST not HR?" | PH violated; HR is misleading time-average; RMST is hazard-free in interpretable time units (Royston-Parmar 2013) |
| "How was tau chosen?" | Pre-specified in SAP based on clinical horizon AND statistical constraint min(largest follow-up); cite Tian 2020 |
| "Fine-Gray or cause-specific Cox?" | Both reported per Putter 2007: FG for CIF prediction, cause-specific Cox for etiology; cite Andersen-Keiding 2012 |
| "MaxCombo direction validation?" | Pre-specified directional constraint; positive z at late-emphasis weight required; cite Magirr-Burman 2021 |
| "Crossover handling?" | ITT primary (treatment policy); RPSFT or IPCW as sensitivity (hypothetical); cite Robins-Tsiatis |
| "Informative censoring?" | Censoring reasons tabulated by arm; symmetric -> no concern; differential -> sensitivity under composite or worst-case |
| "ADTTE CNSR convention?" | Explicit conversion documented: ADaM CNSR=0 means event; converted to event=1 for downstream R/Python |
| "Schoenfeld SS under expected NPH?" | Simulation-based SS using expected hazard pattern (Lakatos 1988); NOT Schoenfeld formula |
| "Multi-state vs Fine-Gray?" | Multi-state when alive->relapse->death is the disease model; Fine-Gray only for CIF prediction of one cause |

## References

- Andersen PK, Keiding N. 2012. Interpretability and importance of functionals in competing risks and multistate models. *Stat Med* 31:1074-1088.
- Andersen PK, Klein JP. 2007. Regression analysis for the cumulative incidence functions. *Scand J Stat* 34:3-16.
- Buyse M, Molenberghs G. 1998. Criteria for the validation of surrogate endpoints in randomized experiments. *Biometrics* 54:1014-1029.
- Cox DR. 1972. Regression models and life-tables. *JRSS-B* 34:187-220.
- Fine JP, Gray RJ. 1999. A proportional hazards model for the subdistribution of a competing risk. *JASA* 94:496-509.
- Karrison TG. 2016. Versatile tests for comparing survival curves based on weighted log-rank statistics. *Stat J* 16:678-690.
- Lakatos E. 1988. Sample sizes based on the log-rank statistic in complex clinical trials. *Biometrics* 44:229-241.
- Lewis EF et al. 2023. Time-to-event endpoints in oncology under the ICH E9(R1) estimand framework. *Pharm Stat* (specific cite varies; see oncology estimand working-group publications).
- Magirr D, Burman CF. 2021. The strong null hypothesis and the MaxCombo test. *Stat Biopharm Res* 15(2):293-298. (Earlier "Cherry-picking in survival analysis" attribution was a blog post, Oct 2022.)
- Mantel N. 1966. Evaluation of survival data and two new rank order statistics arising in its consideration. *Cancer Chemotherapy Reports* 50:163-170.
- Peduzzi P et al. 1995. Importance of events per independent variable in proportional hazards analysis. *J Clin Epidemiol* 48:1503-1510.
- Putter H, Fiocco M, Geskus RB. 2007. Tutorial in biostatistics: competing risks and multi-state models. *Stat Med* 26:2389-2430.
- Putter H, Schumacher M, van Houwelingen HC. 2020. On the relation between the cause-specific hazard and the subdistribution rate for competing risks data: the Fine-Gray model revisited. *Biom J* 62:790-807.
- Royston P, Parmar MKB. 2013. Restricted mean survival time: an alternative to the hazard ratio. *BMC Med Res Methodol* 13:152.
- Schoenfeld DA. 1981. The asymptotic properties of nonparametric tests for comparing survival distributions. *Biometrika* 68:316-319.
- Sun J. 2006. *The Statistical Analysis of Interval-Censored Failure Time Data*. Springer.
- Grambsch PM, Therneau TM. 1994. Proportional hazards tests and diagnostics based on weighted residuals. *Biometrika* 81:515-526.
- Tian L et al. 2020. Empirical comparison of the restricted mean survival time. *Biostatistics*.
- Uno H et al. 2014. Moving beyond the hazard ratio in quantifying the between-group difference in survival analysis. *JCO* 32:2380-2385.
- Xu R, O'Quigley J. 2000. Estimating average regression effect under non-proportional hazards. *Biostatistics* 1:423-439.

## Related Skills

- clinical-biostatistics/effect-measures - HR vs RMST as effect measures; CI methods
- clinical-biostatistics/trial-reporting - ICH E9(R1) estimand framework for time-to-event
- clinical-biostatistics/missing-data-sensitivity - Informative censoring and tipping-point for survival
- clinical-biostatistics/cdisc-data-handling - ADTTE structure and CNSR convention
- clinical-biostatistics/subgroup-analysis - Subgroup HTE for survival endpoints
- clinical-biostatistics/power-and-sample-size - Schoenfeld and Lakatos sample size for TTE
- clinical-biostatistics/multiplicity-graphical - Co-primary survival endpoints
- machine-learning/survival-analysis - Predictive survival models and ML extensions
