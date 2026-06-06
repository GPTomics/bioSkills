---
name: bio-machine-learning-omics-classifiers
description: Builds diagnostic and prognostic classifiers on omics feature matrices with regularized logistic regression, random forest, and gradient-boosted trees, handling the p>>n regime, batch shortcut learning, class imbalance, and probability calibration. Use when building a classifier from expression, methylation, or variant data, choosing an algorithm for high-dimensional small-n data, or diagnosing a suspiciously perfect AUC. For unbiased evaluation see machine-learning/model-validation; for feature selection see machine-learning/biomarker-discovery; for time-to-event outcomes see machine-learning/survival-analysis.
tool_type: python
primary_tool: sklearn
---

## Version Compatibility

Reference examples tested with: pandas 2.2+, scikit-learn 1.4+, xgboost 2.0+, imbalanced-learn 0.12+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

Two high-risk drifts: XGBoost moved `early_stopping_rounds` from `fit()` to the constructor (deprecated 1.6, removed from `fit()` in 2.1); scikit-learn deprecated `LogisticRegression(penalty=)` in 1.8 (use `l1_ratio`+`C`) and `CalibratedClassifierCV(cv='prefit')` in 1.6 (use `FrozenEstimator`). If code throws TypeError/FutureWarning, switch to the constructor / `l1_ratio` / `FrozenEstimator` form.

# Classification Models for Omics Data

**"Build a classifier from my expression data"** -> Start with a regularized linear model (often the ceiling in p>>n), check for batch shortcuts, and treat the probability -- not the label -- as the product.
- Linear (often best): `LogisticRegression(penalty='elasticnet', solver='saga')`
- Trees when nonlinear/interaction signal: `RandomForestClassifier`, `xgboost.XGBClassifier`
- Imbalance: `class_weight='balanced'` or threshold tuning, NOT SMOTE for risk models

## The Single Most Important Modern Insight -- In p>>n, Simple Often Wins and the Probability, Not the Label, Is the Product

Omics classification almost always lives in p>>n (thousands of features, tens-to-hundreds of samples). Two counterintuitive consequences follow. First, more flexible is not better: with n in the dozens the variance of a flexible learner dominates, the full covariance is singular so QDA/full-LDA are undefined, and simple diagonal/linear methods match or beat elaborate ones (Dudoit 2002). "Random forest is the obvious choice for expression data" is a myth -- SVM/regularized logistic frequently win on microarray-style problems (Statnikov 2008), and gradient-boosted trees beat deep nets on tabular/omics data (Grinsztajn 2022). Regularization is the load-bearing wall, not a tuning nicety.

Second, in diagnostic/prognostic use the probability is the product, not the label -- which makes calibration, not accuracy, the thing that breaks silently (Van Calster 2019). A model can rank perfectly (AUC 0.9) and still output dishonest risks. And the most common cause of a beautiful AUC is not skill but a batch artifact: if batch correlates with the outcome, the classifier learns the cleaner technical signal and the performance collapses on any independent cohort.

## Algorithm Choice for p>>n

| Model | Wins when | Overfits / fails when | Calibration | Scaling |
|-------|-----------|------------------------|-------------|---------|
| L1 logistic (lasso) | Sparse signal, want a small signature | Correlated features -> unstable selection; >n true signals | Good (proper loss); shrinks toward base rate | Standardize |
| L2 / elastic-net logistic | Many small correlated effects; omics default | Needs C (and l1_ratio) tuning | Good; preferred when calibration matters | Standardize |
| DLDA / nearest-centroid | Tiny n, roughly linear (Dudoit 2002) | Strong interactions; non-Gaussian | Crude; recalibrate | Variance-scaled |
| Linear SVM | High-dim linear separability (Statnikov 2008) | Heavy overlap; needs C | No native probabilities -- Platt-scale `decision_function` | Critical |
| Random forest | Nonlinear/interaction signal; robust baseline | Sparse-linear signal; tiny n; OOB-as-test leakage | Bagged votes bounded away from 0 and 1 | Scale-invariant |
| GBDT (XGBoost/LightGBM) | Best general tabular performer | Tiny n + deep/many rounds; needs early stopping | Log-loss overfitting tends to overconfident extremes | Scale-invariant |
| Tabular deep nets | Very large n; multimodal/transfer | Typical omics n -> loses to GBDT | Variable; often needs temperature scaling | Standardize |

Tree ensembles are often miscalibrated and the direction depends on the learner and loss: classic boosted ensembles push probabilities toward 0.5 (sigmoid distortion; Niculescu-Mizil 2005), bagged forests are comparatively well-calibrated but bound their votes away from 0 and 1, and modern gradient boosting trained to log-loss for many rounds tends to overfit toward overconfident extremes. Check a reliability curve and recalibrate rather than assuming a direction.

## Decision Tree by Scenario

| Scenario | Recommended approach | Why |
|----------|---------------------|-----|
| Default omics classifier, want a signature | Elastic-net logistic | Often the ceiling in p>>n; sparse + grouping; well-calibrated |
| Suspected nonlinear/interaction (epistasis, thresholds) | Random forest then XGBoost with early stopping | Trees capture interactions; benchmark vs the linear baseline |
| Probabilities will drive a clinical decision | Linear model + calibration check; recalibrate if needed | Probability is the product; AUC is blind to calibration |
| Class imbalance | `class_weight='balanced'` or threshold tuning; never SMOTE for risk | Resampling destroys calibration for no AUC gain |
| Mixed continuous + categorical features | `ColumnTransformer` (scale continuous, encode categorical) | Different feature types need different handling |
| Missing values, especially below-detection | XGBoost/LightGBM native NaN handling | Missingness is often informative (MNAR) |
| Considering a deep net | Only at very large n or multimodal/raw inputs | GBDT beats deep on engineered omics matrices (Grinsztajn 2022) |
| Need unbiased performance / nested CV / calibration metrics | -> machine-learning/model-validation | Evaluation is its own discipline |
| Time-to-event outcome | -> machine-learning/survival-analysis | Censoring needs survival models, not classifiers |

## Core Workflow: Regularized Logistic First

**Goal:** A calibrated, interpretable baseline that is often the best omics classifier.

**Approach:** Standardize inside a Pipeline and fit elastic-net logistic with cross-validated penalty; the L2 component keeps correlated genes together, the L1 component yields a sparse signature.

```python
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegressionCV

# saga supports elasticnet; C = 1/lambda (small C = strong shrinkage). Standardize: the penalty is scale-sensitive.
clf = LogisticRegressionCV(penalty='elasticnet', solver='saga', l1_ratios=[0.1, 0.5, 0.9],
                           Cs=20, cv=5, max_iter=10000, class_weight='balanced')
pipe = Pipeline([('scaler', StandardScaler()), ('clf', clf)])
pipe.fit(X_train, y_train)
```

## Tree Ensembles

**Goal:** Capture nonlinear and interaction structure when the linear baseline leaves signal on the table.

**Approach:** Random forest needs no scaling and is a robust baseline; XGBoost needs a low learning rate, shallow depth, and early stopping (set in the constructor in 2.x) to avoid overfitting tiny n.

```python
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier

rf = RandomForestClassifier(n_estimators=500, max_features='sqrt', min_samples_leaf=3,
                            class_weight='balanced', n_jobs=-1, random_state=0)

# XGBoost 2.x: early_stopping_rounds and eval_metric go in the CONSTRUCTOR, not fit().
# scale_pos_weight is omitted on purpose: like resampling, it reweights the prior and
# distorts calibration -- use it only for hard-label problems, not risk models (see Class Imbalance).
xgb = XGBClassifier(n_estimators=2000, learning_rate=0.03, max_depth=4, subsample=0.8,
                    colsample_bytree=0.5, reg_lambda=1.0,
                    early_stopping_rounds=50, eval_metric='aucpr', n_jobs=-1, random_state=0)
xgb.fit(X_train, y_train, eval_set=[(X_val, y_val)])      # NaN handled natively (missing=np.nan)
```

## Detecting Batch Shortcut Learning

**Goal:** Rule out that a high AUC is a batch artifact rather than biology.

**Approach:** Try to predict the batch from the features and use batch-aware splits; if batch is confounded with the outcome, no correction rescues the design (Soneson 2014) -- fix it at the design stage.

```python
from sklearn.model_selection import cross_val_score, StratifiedGroupKFold
from scipy.stats import chi2_contingency
import pandas as pd

# 1. Can the classifier predict the BATCH? If yes, batch is a strong axis and the label model is suspect.
batch_auc = cross_val_score(pipe, X, batch_labels, cv=5, scoring='roc_auc')
print(f'Batch predictability AUC: {batch_auc.mean():.2f} (high = shortcut risk)')

# 2. Is the outcome associated with batch by design?
print('label vs batch p:', chi2_contingency(pd.crosstab(y, batch_labels))[1])

# 3. Leave-one-batch-out is the honest generalization estimate (usually << random-split CV).
gcv = StratifiedGroupKFold(n_splits=5)
honest = cross_val_score(pipe, X, y, cv=gcv, groups=batch_labels, scoring='roc_auc')
print(f'Batch-aware AUC: {honest.mean():.2f}')
```

## Class Imbalance: What Works and What Fails

**Goal:** Handle a rare positive class without destroying the probabilities.

**Approach:** For a risk model, do not resample -- class-weight cautiously or tune the threshold on a validation fold; SMOTE/oversampling change the training prior, inflate minority probabilities, give no AUC gain, and the same sensitivity is recoverable by moving the threshold (van den Goorbergh 2022). When resampling is unavoidable (a hard-label problem), use an `imblearn` Pipeline so only training folds are resampled.

```python
from imblearn.pipeline import Pipeline as ImbPipeline   # NOT sklearn's Pipeline
from imblearn.over_sampling import SMOTE
from sklearn.linear_model import LogisticRegression

# Correct placement: SMOTE's fit_resample runs only during fit on the train fold, no-op on transform.
imb = ImbPipeline([('smote', SMOTE(random_state=0)),
                   ('clf', LogisticRegression(max_iter=5000))])
# Prefer for risk models: no resampling, then pick the operating threshold by cost on a validation fold.
```

## Probability Calibration

**Goal:** Ensure a "0.9" means a 90% risk, not just a high rank.

**Approach:** Tree ensembles are often miscalibrated -- bagged forests bound votes away from 0 and 1, classic boosting is sigmoid-distorted toward 0.5 (Niculescu-Mizil 2005), and log-loss GBDT can overfit to overconfident extremes -- so check a reliability curve and recalibrate on a disjoint fold. Logistic regression optimizes a proper scoring rule and is usually best-calibrated out of the box. See machine-learning/model-validation for reliability curves, Brier, and the full protocol.

```python
from sklearn.calibration import CalibratedClassifierCV
from sklearn.frozen import FrozenEstimator               # sklearn >=1.6; cv='prefit' deprecated
calibrated = CalibratedClassifierCV(FrozenEstimator(rf.fit(X_tr, y_tr)), method='isotonic')
calibrated.fit(X_cal, y_cal)                             # X_cal disjoint from train and test
```

## Preprocessing, Scaling, Encoding, Missing Data

- Inside the CV fold (refit per fold via Pipeline): feature selection, scaling, log/VST/quantile normalization, ComBat batch correction, imputation, PCA. Fitting any of these on the full matrix is leakage (machine-learning/model-validation).
- Scale-sensitive: regularized logistic, SVM, kNN, PCA, neural nets -- standardize. Scale-invariant: trees/RF/GBDT.
- Variant features: additive 0/1/2 ordinal for trees/additive logistic; one-hot for non-additive (dominant/recessive) effects; CatBoost-style encoding for high-cardinality (HLA).
- Missing data: XGBoost/LightGBM learn a default split direction for NaN; in omics missingness is often MNAR (below detection limit) and a missingness indicator can carry signal -- naive zero-fill conflates "absent" with "not measured."

## Hyperparameters That Matter

| Model | Tune these | Leave default |
|-------|-----------|---------------|
| Logistic | `C` (log-spaced), `l1_ratio`, `class_weight` | solver (saga for elasticnet) |
| Random forest | `max_features`, `min_samples_leaf`, `max_depth` (cap for tiny n) | `n_estimators` (more is safe; 500-1000) |
| XGBoost | `learning_rate`+`n_estimators`+early stopping, `max_depth` (3-6), `subsample`, `colsample_bytree`, `reg_lambda` | most others |

## Per-Method Failure Modes

### Beautiful AUC that is a batch artifact
- **Trigger:** Cases and controls processed in different batches/sites/times.
- **Mechanism:** The classifier exploits the cleaner technical signal; even nested CV is optimistic, and ComBat cannot rescue a confounded design (Soneson 2014).
- **Symptom:** Near-perfect CV AUC; collapse on an independent cohort; the model predicts batch easily.
- **Fix:** Detect with the batch-prediction check; use leave-one-batch-out; fix at design (balance batches across outcome).

### RF/boosting probabilities trusted as risks
- **Trigger:** Reading `predict_proba` from RF or XGBoost as a calibrated risk.
- **Mechanism:** RF bounds votes away from 0 and 1; classic boosting is sigmoid-distorted (Niculescu-Mizil 2005) and log-loss GBDT can overfit to overconfident extremes.
- **Symptom:** Good AUC, reliability curve far from diagonal.
- **Fix:** Recalibrate on a disjoint fold; or prefer logistic when the probability matters.

### SMOTE-before-split / SMOTE for a risk model
- **Trigger:** Resampling before the CV split, or to "fix" imbalance for a probability model.
- **Mechanism:** Synthetic points derived from test samples leak; resampling inflates minority risk and wrecks calibration (van den Goorbergh 2022; Carriero 2025).
- **Symptom:** Inflated CV performance; predicted risks systematically too high.
- **Fix:** `imblearn` Pipeline (train-fold only); for risk models, do not resample -- tune the threshold.

### OOB error read as an unbiased test estimate
- **Trigger:** Reporting RF out-of-bag error after selecting features on the full data.
- **Mechanism:** Selection leaked; OOB then reflects the contaminated feature set.
- **Symptom:** Optimistic OOB; external collapse.
- **Fix:** Selection inside CV; estimate performance by nested CV.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Try a regularized linear model first | Dudoit 2002; Statnikov 2008 | Simple often beats complex in p>>n |
| GBDT over deep nets for tabular omics | Grinsztajn 2022 | Trees handle uninformative features and non-rotational data |
| Do not resample for risk models | van den Goorbergh 2022; Carriero 2025 | Resampling destroys calibration for no AUC gain |
| XGBoost: low LR + many rounds + early stopping | field standard | Prevents overfitting tiny n |
| Report AUPRC + MCC under imbalance | Saito 2015; Chicco 2020 | Accuracy and ROC-AUC mislead when positives are rare |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| XGBoost `early_stopping_rounds` TypeError in `fit()` | Moved to constructor in 2.x | Pass it (and `eval_metric`) in `XGBClassifier(...)` |
| `penalty='l1'` FutureWarning | Deprecated in sklearn 1.8 | Use `l1_ratio=1` + `C` (1.8+) or keep `penalty` on 1.4-1.7 |
| `elasticnet` solver error | Only saga supports it | `solver='saga'` + `l1_ratio` |
| `CalibratedClassifierCV(cv='prefit')` deprecated | sklearn 1.6 | Wrap in `FrozenEstimator` |
| 95% accuracy but useless model | Imbalance + accuracy metric | Report AUPRC/MCC; check the confusion matrix |

## References

- Dudoit S, Fridlyand J, Speed TP. 2002. Comparison of discrimination methods for the classification of tumors using gene expression data. *J Am Stat Assoc* 97:77-87.
- Chawla NV, Bowyer KW, Hall LO, Kegelmeyer WP. 2002. SMOTE: Synthetic Minority Over-sampling Technique. *J Artif Intell Res* 16:321-357.
- Zou H, Hastie T. 2005. Regularization and variable selection via the elastic net. *J R Stat Soc B* 67:301-320.
- Niculescu-Mizil A, Caruana R. 2005. Predicting good probabilities with supervised learning. *Proc 22nd ICML* 625-632.
- Diaz-Uriarte R, Alvarez de Andres S. 2006. Gene selection and classification of microarray data using random forest. *BMC Bioinformatics* 7:3.
- Statnikov A, Wang L, Aliferis CF. 2008. A comprehensive comparison of random forests and support vector machines for microarray-based cancer classification. *BMC Bioinformatics* 9:319.
- Soneson C, Gerster S, Delorenzi M. 2014. Batch effect confounding leads to strong bias in performance estimates obtained by cross-validation. *PLoS ONE* 9:e100335.
- Saito T, Rehmsmeier M. 2015. The precision-recall plot is more informative than the ROC plot when evaluating binary classifiers on imbalanced datasets. *PLoS ONE* 10:e0118432.
- Van Calster B, McLernon DJ, van Smeden M, Wynants L, Steyerberg EW. 2019. Calibration: the Achilles heel of predictive analytics. *BMC Med* 17:230.
- Chicco D, Jurman G. 2020. The advantages of the Matthews correlation coefficient (MCC) over F1 score and accuracy in binary classification evaluation. *BMC Genomics* 21:6.
- Shwartz-Ziv R, Armon A. 2022. Tabular data: deep learning is not all you need. *Inf Fusion* 81:84-90.
- Grinsztajn L, Oyallon E, Varoquaux G. 2022. Why do tree-based models still outperform deep learning on typical tabular data? *NeurIPS Datasets and Benchmarks*.
- van den Goorbergh R, van Smeden M, Timmerman D, Van Calster B. 2022. The harm of class imbalance corrections for risk prediction models. *J Am Med Inform Assoc* 29:1525-1534.
- Carriero A, Luijken K, de Hond A, Moons KGM, van Calster B, van Smeden M. 2025. The harms of class imbalance corrections for machine learning based prediction models. *Stat Med* 44:e10320.

## Related Skills

- machine-learning/model-validation - Nested CV, calibration, and net benefit for the trained classifier
- machine-learning/biomarker-discovery - Select features before modeling (inside the CV fold)
- machine-learning/prediction-explanation - Interpret the classifier and detect shortcuts with SHAP
- machine-learning/survival-analysis - Time-to-event outcomes that classifiers cannot handle
- differential-expression/batch-correction - Batch correction done design-aware, not across the split
- expression-matrix/normalization - Per-sample normalization that is safe outside the CV fold
