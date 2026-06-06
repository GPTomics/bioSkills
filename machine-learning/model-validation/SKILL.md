---
name: bio-machine-learning-model-validation
description: Validates predictive models on omics and biomedical data with nested cross-validation, group/batch/temporal-aware splits, the full data-leakage taxonomy, probability calibration, decision-curve net benefit, optimism correction, sample-size planning, and TRIPOD+AI reporting. Use when estimating model performance honestly, choosing a CV scheme, detecting leakage, or judging whether reported discrimination means the model is actually useful. For feature selection itself see machine-learning/biomarker-discovery; for confirmatory-trial inference see clinical-biostatistics/trial-reporting.
tool_type: python
primary_tool: sklearn
---

## Version Compatibility

Reference examples tested with: numpy 1.26+, scikit-learn 1.4+ (note 1.6/1.8 API changes below).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

scikit-learn drift to watch: `CalibratedClassifierCV(cv='prefit')` was deprecated in 1.6 and removed in 1.8 (it now raises; wrap a fitted model in `sklearn.frozen.FrozenEstimator` instead); `ensemble` default became `'auto'` in 1.6; `method='temperature'` was added in 1.8. If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Model Validation for Biomedical and Omics Data

**"Validate my omics classifier honestly"** -> Keep every data-dependent step inside the resampling loop, never use the same data to both choose and grade, and report calibration and net benefit, not just AUC.
- Nested CV: `GridSearchCV` (inner) wrapped by `cross_val_score` (outer)
- Group/structured: `StratifiedGroupKFold`, `TimeSeriesSplit`
- Calibration: `calibration_curve`, `brier_score_loss`, `CalibratedClassifierCV`

## The Single Most Important Modern Insight -- A Reported Number Is an Honest Estimate Only If Nothing Leaked and Nothing Was Graded on What Was Chosen

A reported performance number is a claim about a data-generating process that will never recur. Almost every inflated result in ML-for-biology traces to one of two root causes: information from the test distribution leaked into model construction, or the same data was used to both choose and grade a decision. A clean train/test split is necessary but nowhere near sufficient -- the leakage has usually already contaminated the test set (a scaler fit on all data, a duplicate patient, ComBat run across the split). Leakage causes a reproducibility crisis across ML-based science (Kapoor 2023), and the bias is largest exactly when the true signal is weakest -- the omics regime.

A second, equally load-bearing insight: discrimination (AUC/C) and calibration (do predicted probabilities match observed frequencies) are *orthogonal*. AUC is invariant to any monotone transform of the score, so it is blind to calibration. For any decision that uses the probability itself, calibration -- not AUC -- is the property that matters, and it is the one routinely ignored (Van Calster 2019, "the Achilles heel").

## Leakage Taxonomy

| Leakage type | How it happens in omics | Symptom | Prevention |
|--------------|-------------------------|---------|------------|
| Preprocessing (most common, most missed) | z-scoring, quantile/library normalization, ComBat/SVA, PCA, kNN/MICE imputation, VST fit on the *full* dataset before splitting | Test performance suspiciously close to train; collapses on external data | Fit every transform inside the CV fold via a `Pipeline` |
| Feature selection (severe special case) | top-k DE genes / highest-variance / univariate filter chosen on all samples, then CV only the classifier | Near-perfect CV from pure noise; unstable selected set | Selection lives in the CV fold (Ambroise 2002) |
| Target / label | a feature is a proxy for or downstream of the outcome (post-diagnosis labs, treatment-derived fields, a collection-site that tracks case/control) | One feature dominates implausibly; fails when removed | Audit temporal/causal admissibility; exclude post-outcome variables |
| Group / patient / replicate | same patient, tumor, organoid, or technical replicate in train and test; `KFold` scatters them | Inflated metrics that vanish under leave-one-group-out | Split by the highest independent unit (`GroupKFold`/`StratifiedGroupKFold`) |
| Batch | batch correlated with outcome and not respected in the split, or ComBat across the train/test boundary | Model discriminates batches not biology; external batch destroys it | Block the split by batch; never run unsupervised correction across the split |
| Temporal | random-splitting time-ordered data; future-period statistics standardize the past | Backtest beats prospective deployment | Time-based split (`TimeSeriesSplit`); never shuffle first |
| Duplicate / homolog | near-identical samples, augmented copies, public-dataset overlap, homologous sequences across the split | Memorization passes as generalization | Deduplicate / cluster-then-split before CV |
| Test-reuse / threshold | repeatedly peeking to pick features, thresholds, "best epoch"; choosing the classification threshold on the test set | Irreproducible SOTA; fragile config | One locked test set; all tuning + thresholds inside nested CV |

## Decision Tree by Scenario

| Scenario / generalization question | Recommended scheme | Why |
|------------------------------------|--------------------|-----|
| "A new sample like training" (and any tuning occurs) | Nested CV: inner `GridSearchCV`, outer `cross_val_score`, Pipeline inside | Tuning and grading on the same CV is optimistic (Cawley-Talbot 2010) |
| "A new patient" (repeated measures) | `GroupKFold`/`StratifiedGroupKFold` by patient/donor | The unit of independence is not the row |
| "A new hospital/site" (transportability) | Leave-one-site-out (internal-external CV) | Approximates external validation |
| "Next year" (time-ordered) | `TimeSeriesSplit` forward-chaining | Random folds leak the future |
| Small n (dozens), need a stable estimate | `RepeatedStratifiedKFold` (5x10) with an interval | A single CV is one high-variance draw |
| Probabilities will drive a decision | Add calibration + decision-curve net benefit | AUC is blind to calibration and utility |
| Final evidence for a clinical model | External/temporal validation + TRIPOD+AI report | Internal CV cannot detect a whole-dataset confound |
| Choosing the features themselves | -> machine-learning/biomarker-discovery | Selection is its own discipline (run it inside the fold) |
| Confirmatory trial inference (HR, p-value) | -> clinical-biostatistics/trial-reporting | Estimand is a treatment effect, not a prediction |

## Nested Cross-Validation

**Goal:** Estimate the performance of the whole procedure (tuning + fit) without optimistic bias.

**Approach:** The inner loop does all tuning, feature selection, and threshold choice; the outer loop grades the winning configuration once on a fold it never touched. The reported number is the aggregate over outer folds; it answers "if I run this pipeline on new data, what do I get?"

```python
from sklearn.model_selection import cross_val_score, StratifiedKFold, GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression

pipe = Pipeline([('scaler', StandardScaler()),
                 ('select', SelectKBest(f_classif)),         # re-fit per inner fold -> no leakage
                 ('clf', LogisticRegression(max_iter=5000))])
grid = {'select__k': [10, 50, 200], 'clf__C': [0.01, 0.1, 1]}

inner = StratifiedKFold(5, shuffle=True, random_state=0)
outer = StratifiedKFold(5, shuffle=True, random_state=1)
search = GridSearchCV(pipe, grid, cv=inner, scoring='roc_auc')
scores = cross_val_score(search, X, y, cv=outer, scoring='roc_auc')   # unbiased estimate
print(f'Nested AUC: {scores.mean():.3f} +/- {scores.std():.3f}')
```

Nested CV is needed whenever model selection happens -- even informal "tried three options, kept the best." Flat CV with tuning is a known reviewer red flag (Varma-Simon 2006).

## Group-Aware, Structured, and Small-Sample CV

**Goal:** Match the CV scheme to the real unit of independence and get a variance-aware estimate.

**Approach:** Pass a grouping vector so no group spans folds; for tiny n, repeat stratified k-fold and report the spread, not a bare number. Standard `KFold` assumes i.i.d. rows, which biomedical data almost never satisfy.

```python
from sklearn.model_selection import StratifiedGroupKFold, RepeatedStratifiedKFold, cross_val_score

groups = meta['patient_id'].values                          # multiple samples per patient
gcv = StratifiedGroupKFold(n_splits=5)                      # group-disjoint AND class-balanced
g_auc = cross_val_score(pipe, X, y, cv=gcv, groups=groups, scoring='roc_auc')

rcv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=0)
r_auc = cross_val_score(pipe, X, y, cv=rcv, scoring='roc_auc')   # report an interval
```

Leave-one-out is high-variance and degenerate for ranking metrics (AUC is undefined on a size-1 test fold) -- prefer repeated stratified k-fold. The .632+ bootstrap (Efron-Tibshirani 1997) is a defensible alternative but is optimistic for zero-apparent-error learners; for internal validation of a single fixed model, bootstrap optimism correction is the cleaner choice.

## Calibration vs Discrimination

**Goal:** Verify that predicted probabilities mean what they say, not just that they rank correctly.

**Approach:** Plot a reliability curve, score it with the proper Brier score, and recalibrate on a held-out fold if needed. AUC measures only ranking; the calibration slope (<1 signals overfitting) and the reliability curve localize the failure.

```python
from sklearn.calibration import calibration_curve, CalibratedClassifierCV
from sklearn.metrics import brier_score_loss
from sklearn.frozen import FrozenEstimator                  # sklearn >=1.6

prob_true, prob_pred = calibration_curve(y_test, p_test, n_bins=10, strategy='quantile')
brier = brier_score_loss(y_test, p_test)                    # proper score: calibration + refinement

# Recalibrate a fitted model on a disjoint calibration fold (cv='prefit' deprecated in 1.6, removed in 1.8):
calibrated = CalibratedClassifierCV(FrozenEstimator(fitted_model), method='isotonic')
calibrated.fit(X_cal, y_cal)                                # X_cal disjoint from train and test
```

Calibration cautions: use `strategy='quantile'` (equal-mass bins) under imbalance; do not report a single Expected Calibration Error as ground truth -- equal-width ECE is biased and reports error even for perfectly calibrated models (Roelofs 2022). Use Platt (`method='sigmoid'`) for small calibration sets, isotonic for hundreds-plus points. Recalibrating on the test set is leakage.

**Net benefit / Decision Curve Analysis (Vickers-Elkin 2006):** `net_benefit = TP/n - (FP/n)*(pt/(1-pt))`, where the threshold probability `pt` encodes the relative harm of a false positive. Plot it against treat-all and treat-none references; a model is clinically useful only where it sits above both. DCA requires good calibration to be valid and is the bridge from statistical performance to clinical usefulness -- a model can have high AUC yet zero net benefit at every plausible threshold.

## Metric Selection in Imbalanced Data

| Metric | Use | Trap |
|--------|-----|------|
| Accuracy | Almost never headline it under imbalance | At 5% prevalence, "always negative" scores 95% |
| AUC / C | Discrimination, prevalence-independent | Blind to calibration; not a usefulness measure |
| AUPRC (average precision) | Rare-positive problems | Baseline is the prevalence, not 0.5 -- state it (Saito 2015) |
| Brier / log-loss | When probabilities are used | Proper; not comparable across prevalences without scaling |
| MCC | Balanced single-threshold summary | Still threshold-dependent (Chicco 2020) |
| F1 | Retrieval-style problems | Ignores true negatives; assumes a cost ratio |

The multiple-threshold problem: reporting the *best* F1/accuracy over thresholds is optimistic, and choosing that threshold on the test set is leakage. Pick the operating point on a separate fold (or by net benefit), then report the locked-threshold metric once; prefer threshold-free curves (ROC, PR, calibration) plus one pre-specified operating point.

## External Validation, Optimism, Sample Size, and TRIPOD+AI

- **Internal vs external.** Internal validation (bootstrap optimism correction, repeated/nested CV) estimates reproducibility on new patients from the *same* source; external validation (different time, place, setting) estimates transportability and is the usual point of failure -- calibration degrades first (slope <1, intercept shift). Internal-external CV (leave-one-cluster-out) is the recommendation when multiple cohorts exist (Steyerberg 2001).
- **Optimism and shrinkage.** Apparent performance overstates the future; the gap (optimism) grows with more predictors, more flexibility, smaller n. Remedy with a uniform shrinkage factor (the bootstrap calibration slope) or penalized estimation. A development-data calibration slope <1 *is* the optimism signal.
- **Sample size.** The "10 events per variable" heuristic (Peduzzi 1996) is obsolete; the standard is Riley et al.'s minimum-sample-size framework (2019, Stat Med Parts I-II), which sizes for shrinkage >=0.9 and precise risk estimation (`pmsampsize`). For p>>n omics these formulas are out of regime, which is precisely why heavy penalization + nested validation, not unpenalized multivariable fits, are mandatory.
- **Reporting.** TRIPOD+AI (Collins 2024, *BMJ* 385:e078378) supersedes TRIPOD 2015 and is the 2024+ target for any biomedical predictive-model claim -- it demands data-splitting and leakage controls, calibration (not just discrimination), fairness/subgroup performance, and uncertainty. PROBAST+AI is the companion risk-of-bias appraisal.

## Per-Method Failure Modes

### Preprocessing fit before the split
- **Trigger:** `StandardScaler().fit_transform(X)` (or ComBat, PCA, imputation) on all data, then CV.
- **Mechanism:** The fitted parameters encode the test rows.
- **Symptom:** Test variance tiny; drop on external data.
- **Fix:** Put every transform in the Pipeline so `fit` only sees training folds.

### Threshold or best-of-many chosen on the test set
- **Trigger:** Reporting the best F1 over thresholds, or the best of several CV runs.
- **Mechanism:** Each peek leaks; over many tries the test set becomes a training set.
- **Symptom:** Irreproducible "SOTA"; a fresh test set disappoints.
- **Fix:** Lock one test set, pre-specify metric and threshold rule, choose thresholds on a separate fold.

### SMOTE/resampling to fix imbalance breaks calibration
- **Trigger:** Oversampling/SMOTE for a *risk* model.
- **Mechanism:** Changing training prevalence inflates minority-class probabilities; no AUC gain (van den Goorbergh 2022).
- **Symptom:** Good AUC, badly miscalibrated risks.
- **Fix:** Do not resample for probability models; move the threshold on a calibrated model. If resampled, use `imblearn.pipeline.Pipeline` (train-fold only).

### LOO for a ranking metric
- **Trigger:** Leave-one-out with AUC.
- **Mechanism:** AUC is undefined within a size-1 fold; pooling OOF predictions then scoring once is not equivalent to averaging fold scores for non-decomposable metrics.
- **Symptom:** Unstable or misleading AUC.
- **Fix:** Use repeated stratified k-fold; reserve `cross_val_predict` for visuals, not the headline metric.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Selection/preprocessing inside every fold; nested CV for tuning | Ambroise 2002; Varma-Simon 2006 | Same-data tune-and-grade is optimistic |
| Repeated 5-fold x ~10, report the spread | field standard | A single CV is one high-variance draw at small n |
| Calibration slope ~1; <1 means overfitting | Van Calster 2019 | Basis of shrinkage |
| Sample size from Riley framework (shrinkage >=0.9) | Riley 2019 | "10 EPV" is obsolete |
| Report per TRIPOD+AI | Collins 2024 | 2024+ standard: discrimination + calibration + fairness |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `cv='prefit'` warns or errors | Deprecated in 1.6, removed in 1.8 (now raises) | Wrap the fitted model in `FrozenEstimator` |
| Calibration looks perfect on the test set | Calibrated on the evaluation data | Calibrate on a disjoint fold |
| Scaler/selector fit outside the Pipeline | Preprocessing leakage | Move into the Pipeline passed to `cross_val_score`/`GridSearchCV` |
| `groups=` ignored | Not threaded to the splitter | Pass `groups=` to `cross_validate`/`GridSearchCV.fit` |
| `cross_val_predict` used as the headline AUC | Non-decomposable metric over pooled OOF | Average per-fold scores instead |

## References

- Peduzzi P, Concato J, Kemper E, Holford TR, Feinstein AR. 1996. A simulation study of the number of events per variable in logistic regression analysis. *J Clin Epidemiol* 49:1373-1379.
- Efron B, Tibshirani R. 1997. Improvements on cross-validation: the .632+ bootstrap method. *J Am Stat Assoc* 92:548-560.
- Steyerberg EW, Harrell FE, Borsboom GJ, et al. 2001. Internal validation of predictive models. *J Clin Epidemiol* 54:774-781.
- Ambroise C, McLachlan GJ. 2002. Selection bias in gene extraction on the basis of microarray gene-expression data. *PNAS* 99:6562-6566.
- Vickers AJ, Elkin EB. 2006. Decision curve analysis: a novel method for evaluating prediction models. *Med Decis Making* 26:565-574.
- Varma S, Simon R. 2006. Bias in error estimation when using cross-validation for model selection. *BMC Bioinformatics* 7:91.
- Cawley GC, Talbot NLC. 2010. On over-fitting in model selection and subsequent selection bias in performance evaluation. *J Mach Learn Res* 11:2079-2107.
- Saito T, Rehmsmeier M. 2015. The precision-recall plot is more informative than the ROC plot when evaluating binary classifiers on imbalanced datasets. *PLoS ONE* 10:e0118432.
- Van Calster B, McLernon DJ, van Smeden M, Wynants L, Steyerberg EW. 2019. Calibration: the Achilles heel of predictive analytics. *BMC Med* 17:230.
- Riley RD, Snell KIE, Ensor J, et al. 2019. Minimum sample size for developing a multivariable prediction model: Parts I-II. *Stat Med* 38:1262-1296.
- Chicco D, Jurman G. 2020. The advantages of the Matthews correlation coefficient (MCC) over F1 score and accuracy in binary classification evaluation. *BMC Genomics* 21:6.
- Roelofs R, Cain N, Shlens J, Mozer MC. 2022. Mitigating bias in calibration error estimation. *Proc AISTATS* PMLR 151:4036-4054.
- van den Goorbergh R, van Smeden M, Timmerman D, Van Calster B. 2022. The harm of class imbalance corrections for risk prediction models. *J Am Med Inform Assoc* 29:1525-1534.
- Whalen S, Schreiber J, Noble WS, Pollard KS. 2022. Navigating the pitfalls of applying machine learning in genomics. *Nat Rev Genet* 23:169-181.
- Kapoor S, Narayanan A. 2023. Leakage and the reproducibility crisis in machine-learning-based science. *Patterns* 4:100804.
- Collins GS, Moons KGM, Dhiman P, et al. 2024. TRIPOD+AI statement. *BMJ* 385:e078378.

## Related Skills

- machine-learning/biomarker-discovery - Feature selection run inside the CV fold
- machine-learning/omics-classifiers - Model training, calibration directions, and imbalance handling
- machine-learning/survival-analysis - Validation metrics for time-to-event models
- experimental-design/batch-design - Designing out batch-outcome confounding before analysis
- experimental-design/multiple-testing - FDR control for high-dimensional testing
- clinical-biostatistics/trial-reporting - Confirmatory-trial reporting and the prediction-vs-inference boundary
