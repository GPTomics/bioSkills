---
name: bio-machine-learning-survival-analysis
description: Builds and validates predictive time-to-event models on clinical and omics data with penalized Cox, random survival forests, gradient-boosted and deep survival models, and prediction-grade evaluation (Uno's C, time-dependent AUC, integrated Brier, calibration, competing risks). Use when building an individualized risk predictor or prognostic omics signature, choosing a survival model, or evaluating one beyond the C-index. For Kaplan-Meier, log-rank, and classical Cox hazard-ratio inference in a trial see clinical-biostatistics/survival-analysis.
tool_type: python
primary_tool: scikit-survival
---

## Version Compatibility

Reference examples tested with: scikit-survival 0.22+, lifelines 0.30+, numpy 1.26+, pandas 2.2+ (pycox 0.3+ for deep models).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

scikit-survival requires the target `y` to be a structured array with a boolean event field and a float time field; its IPCW metrics need the *training* `y` first. lifelines `concordance_index` expects higher-score = longer-survival (negate the partial hazard). If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Predictive Survival Modeling

**"Build a validated risk model from time-to-event data"** -> Fit a penalized Cox or ensemble survival model, then evaluate with censoring-robust discrimination AND calibration, not the C-index alone.
- Penalized Cox / RSF / boosting: `sksurv.linear_model.CoxnetSurvivalAnalysis`, `sksurv.ensemble.RandomSurvivalForest`
- Evaluation: `concordance_index_ipcw`, `cumulative_dynamic_auc`, `integrated_brier_score`
- Deep survival (large n): `pycox` DeepSurv / DeepHit

## The Single Most Important Modern Insight -- C-Index Only Is the Cardinal Sin

The C-index is necessary but radically insufficient, for three reasons most papers miss. It is **censoring-distribution-dependent**: Harrell's C is biased upward under heavy censoring and gives different values on cohorts that differ only in follow-up -- report Uno's IPCW C with an explicit truncation tau instead (Uno 2011). It is **invariant to any monotone transform of the risk score**, so a model can have an excellent C and be wildly miscalibrated and clinically harmful -- C measures ranking, not the correctness of the predicted probabilities. And it is **insensitive**: adding a genuinely useful marker barely moves it, which is exactly why reclassification metrics (NRI/IDI) were invented. Decision-grade evaluation is Uno's C(tau) + time-dependent AUC(t) + integrated Brier vs a Kaplan-Meier baseline + calibration curves, all on honestly held-out or external data.

## ML-vs-Confirmatory Scope Boundary

Two survival cultures, two skills; mixing them is the most common authoring error. Resolve every overlap by one question: *is the goal to estimate and test a treatment effect in a (pre-specified) study, or to build and validate a risk-prediction model?*

| This skill (machine-learning) -- prediction | clinical-biostatistics/survival-analysis -- inference |
|---------------------------------------------|--------------------------------------------------------|
| Estimand is an individualized risk (survival curve, risk score, CIF) | Estimand is a treatment effect (hazard ratio with CI, p-value) |
| Penalized Cox, RSF, boosting, DeepSurv/DeepHit | Kaplan-Meier, log-rank, classical low-dimensional Cox |
| p>>n omics signatures; feature selection inside the resampling loop | Pre-specified analysis plan, FWER control, regulated trial |
| Judged by out-of-sample prediction (Uno's C, IBS, calibration) | Judged by validity of the HR and PH diagnostics (`cox.zph`) |
| PH is an assumption to relax (RSF/DeepHit do not assume it) | PH is a hypothesis whose violation invalidates the reported HR |

PH concepts, censoring definitions, and the Cox partial likelihood are foundational to both -- stated in clinical-biostatistics and referenced here. This skill's value begins at "I want a validated predictor."

## Model Taxonomy

| Model | Assumes PH? | Handles p>>n? | Competing risks? | Best when |
|-------|-------------|----------------|-------------------|-----------|
| Penalized Cox (elastic-net, coxnet) | Yes | Yes -- the omics workhorse; elastic-net handles correlated genes | Cause-specific by recoding | Sparse, interpretable, reproducible risk score; the default |
| Random Survival Forest | No | Yes (tune mtry/nodesize) | Yes (per-cause CIF) | Nonlinear/interaction effects, non-PH, moderate n |
| Gradient-boosted survival | Componentwise: yes (sparse); tree base: no | Yes (componentwise selects) | Via cause-specific | Boosting accuracy + sparsity, or relaxing PH with trees |
| Survival SVM | No (optimizes concordance) | Yes (kernel) | No | Pure ranking goals; gives a score, not a survival function |
| DeepSurv (pycox CoxPH) | Yes (NN replaces the linear predictor) | Needs large n | No | Large n, nonlinear main effects, PH plausible |
| DeepHit | No (discrete-time PMF) | Needs large n | Yes -- purpose-built | Large n + competing risks + non-PH |
| Cox-Time | No (time-dependent NN) | Needs large n | Discrete-hazard extensions | Large n, non-PH, flexible survival function |

Load-bearing empirical fact: on typical clinical-omics n, **penalized Cox and RSF are very hard to beat**; deep survival models usually only pull ahead at large n and/or with genuine non-PH or competing-risk structure. Start with elastic-net Cox and RSF baselines; escalate to deep models only if they demonstrably beat those on a held-out set.

## Decision Tree by Scenario

| Scenario | Recommended approach | Why |
|----------|---------------------|-----|
| Prognostic omics signature, p>>n | Elastic-net Cox (coxnet), selection inside nested CV | Sparse + correlated-gene grouping; the workhorse |
| Nonlinear/interaction effects, non-PH suspected | Random survival forest | Assumes no PH; captures interactions |
| Large n, suspected nonlinear main effects | DeepSurv, benchmarked against coxnet/RSF | Deep only earns its keep at large n |
| Competing events (death from other causes) | Fine-Gray (CIF) or DeepHit; report cause-specific too | 1-KM overestimates incidence; sHR is not the rate effect |
| Dynamic prediction with updating biomarkers | Landmarking | Internal time-varying covariates cannot be plugged into the future |
| Evaluating any survival model | Uno's C(tau) + AUC(t) + IBS-vs-KM + calibration | C-index alone is insufficient |
| KM curve, log-rank, or a trial hazard ratio | -> clinical-biostatistics/survival-analysis | Confirmatory inference, not prediction |
| Selecting the prognostic genes | -> machine-learning/biomarker-discovery (same irreproducibility) | Selection is its own discipline, inside the loop |

## Fitting Predictive Survival Models

**Goal:** Fit a penalized Cox and an RSF baseline with the correct target format.

**Approach:** Build the structured `y` (boolean event, float time), fit coxnet with an elastic-net mix and a baseline model for survival functions, and an RSF for nonlinear structure.

```python
from sksurv.util import Surv
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.ensemble import RandomSurvivalForest

# event MUST be bool, time float. Field order in from_arrays is (event, time).
y = Surv.from_arrays(event=df['status'].astype(bool), time=df['time'].astype(float))

# Elastic-net Cox: l1_ratio in (0,1]; fit_baseline_model=True enables predict_survival_function.
coxnet = CoxnetSurvivalAnalysis(l1_ratio=0.9, alpha_min_ratio=0.01, fit_baseline_model=True)
coxnet.fit(X_train, y_train)

rsf = RandomSurvivalForest(n_estimators=500, min_samples_leaf=15, max_features='sqrt', n_jobs=-1)
rsf.fit(X_train, y_train)
risk = coxnet.predict(X_test)              # a risk score (higher = higher risk), NOT a probability
```

## Prediction-Grade Evaluation

**Goal:** Report discrimination AND calibration on out-of-sample data, not the C-index alone.

**Approach:** Use Uno's IPCW C (truncated at tau), time-dependent AUC over a horizon grid, integrated Brier vs the KM baseline, and a calibration check at clinical horizons. IPCW metrics take the *training* `y` first to estimate the censoring distribution.

```python
import numpy as np
from sksurv.metrics import concordance_index_ipcw, cumulative_dynamic_auc, integrated_brier_score

c_uno = concordance_index_ipcw(y_train, y_test, risk, tau=t_horizon)[0]      # censoring-robust

times = np.percentile(y_test['time'][y_test['event']], np.linspace(10, 80, 15))   # inside follow-up
auc_t, mean_auc = cumulative_dynamic_auc(y_train, y_test, risk, times)

surv_fns = coxnet.predict_survival_function(X_test)              # needs fit_baseline_model=True
surv_prob = np.vstack([[fn(t) for t in times] for fn in surv_fns])
ibs = integrated_brier_score(y_train, y_test, surv_prob, times)  # compare against the KM-only IBS
print(f"Uno C: {c_uno:.3f}  mean AUC(t): {mean_auc:.3f}  IBS: {ibs:.3f}")
```

Calibration (the most decision-relevant, most-ignored axis): at a clinical horizon, plot predicted P(event by t) against the observed event probability (KM within risk groups, or a smooth calibration curve), and summarize with ICI/E50/E90 (Austin 2020) -- never Hosmer-Lemeshow. Calibration is what usually breaks on external validation even when discrimination is preserved.

## Competing Risks

A competing event precludes the event of interest (death from another cause precludes cause-specific death). Treating competing events as ordinary censoring is wrong: **Kaplan-Meier overestimates cumulative incidence** (1-KM >= the cumulative incidence function), so report the CIF (Aalen-Johansen), not 1-KM. Two hazards answer two questions: the **cause-specific hazard** (Cox censoring competing events) is the etiologic rate; the **Fine-Gray subdistribution hazard** maps to the CIF and is the prognostic/absolute-risk target -- but a Fine-Gray sHR is NOT the cause-specific HR and not an effect on the event rate (a covariate can raise a CIF purely by lowering the competing hazard). Since this skill is the prediction regime, the CIF is usually the target; report both models for a complete picture. DeepHit is purpose-built for competing risks at large n; evaluate with competing-risks-aware concordance (Wolbers 2009) and CIF-based Brier.

## Censoring, Immortal Time, and Landmarking

- **Non-informative censoring** is assumed by the Cox likelihood, KM, and all IPCW metrics: censored subjects must be representative of those still at risk. If sicker patients drop out (informative censoring), survival is overestimated and IPCW does not fix it with a misspecified censoring model. Administrative censoring is the benign case; informative censoring is generally untestable and needs sensitivity analysis. Truncate IPCW metrics at tau because tail weights explode.
- **Immortal time bias** -- defining a group by a post-baseline event ("patients who received treatment") manufactures a survival advantage from guaranteed event-free time (Levesque 2010). Endemic in EHR/omics; fix with time-varying exposure, landmarking, or target-trial emulation.
- **Landmarking** (van Houwelingen 2007) -- internal time-varying covariates (a biomarker that changes with disease) cannot be plugged into the future; pick a landmark time, restrict to those still event-free, use covariate values as of the landmark, and predict forward. It is the robust route to honest dynamic prediction.

## High-Dimensional Penalized Cox and Validation

In p>>n, elastic-net Cox beats LASSO for stability with correlated genes (LASSO arbitrarily keeps one of a correlated group). The same irreproducibility as biomarker discovery applies: different cohorts select near-disjoint gene sets at similar performance, so the deliverable is the *prediction*, not the gene list. Feature selection and tuning MUST live inside the resampling loop -- nested CV with survival metrics (Uno's C, IBS) as the objective; selecting genes on the full data then CV-ing the final model is leakage producing grossly optimistic signatures. Optimism-correct internal validation via Harrell's bootstrap, validate externally, and report per TRIPOD+AI (Collins 2024).

## Per-Method Failure Modes

### Reporting only the C-index
- **Trigger:** Summarizing a survival model by Harrell's C alone.
- **Mechanism:** C is censoring-dependent, monotone-invariant (blind to calibration), and insensitive.
- **Symptom:** Great C, badly miscalibrated absolute risks; ranking-equivalent models look identical.
- **Fix:** Uno's C(tau) + AUC(t) + IBS-vs-KM + calibration curves on out-of-sample data.

### Kaplan-Meier under competing risks
- **Trigger:** Using 1-KM for incidence when a competing event exists.
- **Mechanism:** 1-KM assumes competing-event subjects could still have the event; they cannot.
- **Symptom:** Incidence overestimated; sums across causes exceed 1.
- **Fix:** Report the CIF (Aalen-Johansen); use cause-specific and Fine-Gray models.

### Misverbalizing a Fine-Gray coefficient
- **Trigger:** Saying a Fine-Gray sHR "increases the rate of the event."
- **Mechanism:** The subdistribution hazard maps to cumulative incidence, not the event rate.
- **Symptom:** Causal/etiologic claims from a prognostic model.
- **Fix:** State it as an effect on *cumulative incidence*; report cause-specific too.

### Immortal time bias
- **Trigger:** Grouping subjects by a post-baseline event.
- **Mechanism:** Guaranteed event-free time is attributed to the exposed group.
- **Symptom:** Spectacular development performance, external collapse.
- **Fix:** Time-varying exposure, landmarking, or target-trial emulation.

### Selection-before-CV in a signature
- **Trigger:** Selecting genes on all data, then CV-ing the final Cox model.
- **Mechanism:** Held-out folds informed selection (the dominant overfitting capacity in p>>n).
- **Symptom:** Optimistic, irreproducible signature.
- **Fix:** Selection inside nested CV with survival metrics; elastic net for stability.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Report Uno's C with explicit tau | Uno 2011 | Harrell's C is censoring-dependent and biased upward |
| Evaluate beyond C: AUC(t), IBS-vs-KM, calibration | Graf 1999; Austin 2020 | C is monotone-invariant and insensitive |
| CIF (not 1-KM) under competing risks | Putter 2007 | 1-KM overestimates incidence |
| Selection inside nested CV with survival metrics | field standard | Selection-before-CV leaks in p>>n |
| Start with penalized Cox / RSF baselines | neutral benchmarks | Deep models rarely beat them at clinical n |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| sksurv crash / silent misbehavior on `y` | y not a structured (bool event, float time) array | Use `Surv.from_arrays(event=..bool, time=..float)` |
| IPCW metric wrong | Training `y` not passed first | `concordance_index_ipcw(y_train, y_test, risk, tau=)` |
| `predict_survival_function` errors | coxnet without baseline | `CoxnetSurvivalAnalysis(fit_baseline_model=True)` |
| lifelines C reported as 1-C | partial hazard not negated | `concordance_index(time, -predict_partial_hazard(df), event)` |
| pycox survival nonsense | baseline hazard not computed | `model.compute_baseline_hazards()` before `predict_surv_df` |
| `times` for AUC/IBS out of range | beyond largest uncensored test time | Clip `times` to the observed follow-up |

## References

- Harrell FE, Lee KL, Mark DB. 1996. Multivariable prognostic models. *Stat Med* 15:361-387.
- Graf E, Schmoor C, Sauerbrei W, Schumacher M. 1999. Assessment and comparison of prognostic classification schemes for survival data. *Stat Med* 18:2529-2545.
- Fine JP, Gray RJ. 1999. A proportional hazards model for the subdistribution of a competing risk. *J Am Stat Assoc* 94:496-509.
- Heagerty PJ, Zheng Y. 2005. Survival model predictive accuracy and ROC curves. *Biometrics* 61:92-105.
- Putter H, Fiocco M, Geskus RB. 2007. Tutorial in biostatistics: competing risks and multi-state models. *Stat Med* 26:2389-2430.
- van Houwelingen HC. 2007. Dynamic prediction by landmarking in event history analysis. *Scand J Stat* 34:70-85.
- Ishwaran H, Kogalur UB, Blackstone EH, Lauer MS. 2008. Random survival forests. *Ann Appl Stat* 2:841-860.
- Wolbers M, Koller MT, Witteman JCM, Steyerberg EW. 2009. Prognostic models with competing risks. *Epidemiology* 20:555-561.
- Levesque LE, Hanley JA, Kezouh A, Suissa S. 2010. Problem of immortal time bias in cohort studies. *BMJ* 340:b5087.
- Simon N, Friedman J, Hastie T, Tibshirani R. 2011. Regularization paths for Cox's proportional hazards model via coordinate descent. *J Stat Softw* 39:1-13.
- Uno H, Cai T, Pencina MJ, D'Agostino RB, Wei LJ. 2011. On the C-statistics for evaluating overall adequacy of risk prediction procedures with censored survival data. *Stat Med* 30:1105-1117.
- Katzman JL, Shaham U, Cloninger A, et al. 2018. DeepSurv: personalized treatment recommender system using a Cox proportional hazards deep neural network. *BMC Med Res Methodol* 18:24.
- Lee C, Zame WR, Yoon J, van der Schaar M. 2018. DeepHit: a deep learning approach to survival analysis with competing risks. *Proc AAAI* 32:2314-2321.
- Austin PC, Harrell FE, van Klaveren D. 2020. Graphical calibration curves and the integrated calibration index (ICI) for survival models. *Stat Med* 39:2714-2742.
- Polsterl S. 2020. scikit-survival: a library for time-to-event analysis built on top of scikit-learn. *J Mach Learn Res* 21:1-6.
- Collins GS, Moons KGM, Dhiman P, et al. 2024. TRIPOD+AI statement. *BMJ* 385:e078378.

## Related Skills

- clinical-biostatistics/survival-analysis - Kaplan-Meier, log-rank, classical Cox inference, PH diagnostics for trials
- machine-learning/model-validation - Nested CV, calibration, and optimism correction shared with prediction models
- machine-learning/biomarker-discovery - Selecting prognostic genes inside the resampling loop
- differential-expression/de-results - Pre-filter candidate prognostic genes
- clinical-databases/variant-prioritization - Clinical interpretation of prognostic variants
