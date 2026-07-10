---
name: bio-workflows-biomarker-pipeline
description: End-to-end biomarker discovery workflow from expression data to validated biomarker panels. Covers feature selection with Boruta/LASSO, leakage-safe cross-validation, calibration, and SHAP interpretation. Use when building and validating diagnostic or prognostic biomarker signatures from omics data.
tool_type: python
primary_tool: sklearn
workflow: true
depends_on:
  - machine-learning/biomarker-discovery
  - machine-learning/model-validation
  - machine-learning/omics-classifiers
  - machine-learning/prediction-explanation
qc_checkpoints:
  - after_selection: "Selected features 5-200, stability index reported alongside count"
  - after_cv: "Selection inside the CV pipeline; AUC reported with fold spread; AUPRC/MCC if imbalanced"
  - after_interpretation: "SHAP used as a shortcut/batch audit, aggregated over modules, not as the validated panel"
  - after_validation: "Hold-out AUC with bootstrap CI plus calibration (Brier); external cohort for the real bar"
---

## Version Compatibility

Reference examples tested with: numpy 1.26+, pandas 2.2+, scikit-learn 1.4+, shap 0.47+ (the `feature_perturbation='auto'` estimand and per-class 3-D `.values` behavior the code relies on; xgboost 2.0+ optional).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

scikit-learn drift: `CalibratedClassifierCV(cv='prefit')` deprecated in 1.6 (use `FrozenEstimator`); `LogisticRegression(penalty=)` deprecated in 1.8, and `LogisticRegressionCV(penalty='l1')` too -- the 1.8+ migration drops `penalty=` entirely and passes `l1_ratios=(1.0,)` alone (leave `penalty` at its default; `penalty='elasticnet'` still emits the FutureWarning). XGBoost moved `early_stopping_rounds` to the constructor in 2.x. If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Biomarker Discovery Pipeline

**"Build a validated biomarker panel from my omics data"** -> Orchestrate group-aware splitting, feature selection, leakage-safe cross-validation, calibration, and SHAP interpretation to produce a robust, honestly-validated biomarker signature.

This is a workflow skill: it owns the chaining decisions and hand-offs, not the internals of any one step. Every step below cross-references the component skill that teaches its mechanism.

## The governing principle

The whole pipeline stands or falls on four commitments made at the seams; each one, if broken, inflates the reported performance and a held-out set cannot detect the leak because it was already contaminated.

1. **The independent unit of splitting is the highest biological unit — patient/donor/site, NOT the sample — and it is committed first.** Multiple biopsies, longitudinal samples, or technical replicates from one subject in both train and test is group leakage; the model memorizes the subject, not the biology. Split with `GroupKFold`/`StratifiedGroupKFold` on a subject key. For single-cell-derived features the unit is the donor, not the cell.
2. **Every data-dependent transform is fit INSIDE the CV fold** — scaling, library-size/quantile normalization, ComBat/SVA, PCA/UMAP, imputation, AND feature selection. The discovery panel may be selected on all training data (that IS the deliverable), but the performance NUMBER must come from a pipeline that re-runs selection per fold. Selection is the dominant overfitting capacity in p>>n and gives near-perfect apparent accuracy on pure noise (Ambroise & McLachlan 2002).
3. **The locked test set is touched exactly once.** Every threshold, feature count, hyperparameter, and "best epoch" chosen on it leaks; when hyperparameters are tuned, use nested CV to report performance (Varma & Simon 2006).
4. **The metric is matched to the data regime, and calibration is separate from discrimination.** AUC for discrimination, AUPRC/MCC when imbalanced, and Brier + a reliability curve whenever risk estimates will be used — AUC is invariant to any monotone score transform, so it says nothing about calibration.

## Workflow Overview

```
Expression matrix + Metadata
    |
    v
[1. Data Preparation] -----> StandardScaler, train/test split
    |
    v
[2. Feature Selection] ----> Boruta or LASSO stability selection
    |
    v
[3. Model Training] -------> Pipeline with selection inside CV (leakage-safe)
    |
    v
[4. Model Interpretation] -> SHAP values, feature importance
    |
    v
[5. Validation] -----------> Hold-out test, bootstrap CI
    |
    v
Validated biomarker panel + classifier
```

## Step 1: Data Preparation

**Goal:** Load the matrix and hold out a GROUP-aware test set before anything is fit.

**Approach:** Split by the subject key so no subject appears in both train and test, then fit the scaler on training only; per-fold scaling is re-applied inside the CV pipeline in Step 3.

```python
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler

expr = pd.read_csv('expression.csv', index_col=0)
meta = pd.read_csv('metadata.csv', index_col=0)

X = expr.T  # samples x genes
# y must be 0/1: brier_score_loss and calibration_curve raise on string labels, and sklearn orders
# classes alphabetically -- for a case/control column that makes 'control' the positive class, so
# predict_proba[:, 1], the SHAP [:, :, 1] slice, and Brier all silently describe the wrong class.
# AUC is symmetric and will not expose the flip. Encode the disease class as 1 explicitly.
POSITIVE_CLASS = 'disease'
y = (meta.loc[X.index, 'condition'].values == POSITIVE_CLASS).astype(int)
# The critical key: the SUBJECT (patient/donor/site), not the sample. If truly one
# sample per subject, groups = X.index; otherwise it MUST be the subject id.
groups = meta.loc[X.index, 'subject_id'].values

# Group- AND class-aware hold-out: take one StratifiedGroupKFold fold as the test set so no
# subject spans train/test (train_test_split(stratify=y) alone would leak repeated subjects).
sgkf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)   # 1/5 held out (~0.2)
train_idx, test_idx = next(sgkf.split(X, y, groups))
X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
y_train, y_test = y[train_idx], y[test_idx]
groups_train = groups[train_idx]

# Fit scaler on training only to prevent data leakage
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)
```

**QC Checkpoint 1:** Check class balance, sample counts, and group separation
- Minimum 10 samples per class recommended; classes reasonably balanced (ratio <3:1)
- Confirm NO subject id appears in both train and test (`set(groups[train_idx]) & set(groups[test_idx])` is empty)

## Step 2: Feature Selection

**Goal:** Produce the discovery panel (all-relevant with Boruta, or a stable minimal set with LASSO).

**Approach:** Optionally pre-filter, then run the selector and map the mask back to the full feature space for downstream indexing.

### Option A: Boruta (All-Relevant Selection)

```python
import numpy as np
from boruta import BorutaPy
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, f_classif

# Pre-filter if >10k features. selected_idx is a positional boolean mask aligned to X_train.columns.
if X_train_scaled.shape[1] > 10000:
    selector = SelectKBest(f_classif, k=5000)
    selector.fit(X_train_scaled, y_train)
    prefilter_idx = np.where(selector.get_support())[0]
    X_train_filt = X_train_scaled[:, prefilter_idx]
else:
    prefilter_idx = None
    X_train_filt = X_train_scaled

# max_depth=5: Shallow trees for stable importances
rf = RandomForestClassifier(n_estimators=100, max_depth=5, n_jobs=-1, random_state=42)
# max_iter=100: Usually sufficient; 200 if many tentative
boruta = BorutaPy(rf, n_estimators='auto', max_iter=100, random_state=42, verbose=0)
boruta.fit(X_train_filt, y_train)

# Map the (possibly pre-filtered) Boruta mask back onto the FULL feature space.
selected_idx = np.zeros(X_train.shape[1], dtype=bool)
selected_idx[prefilter_idx[boruta.support_] if prefilter_idx is not None else boruta.support_] = True
print(f'Selected {selected_idx.sum()} features')
```

### Option B: LASSO Stability Selection

```python
from sklearn.linear_model import LogisticRegressionCV
import numpy as np

# n_bootstrap=100: Quick; use 500 for publication
n_bootstrap = 100
stability_scores = np.zeros(X_train_scaled.shape[1])

for i in range(n_bootstrap):
    idx = np.random.choice(len(y_train), size=len(y_train), replace=True)
    # Cs=10: 10 regularization values to search
    model = LogisticRegressionCV(penalty='l1', solver='saga', Cs=10, cv=3, random_state=i, max_iter=1000)
    model.fit(X_train_scaled[idx], y_train[idx])
    stability_scores += (model.coef_[0] != 0).astype(int)

stability_scores /= n_bootstrap
# stability_threshold=0.6: Standard; 0.8 for strict
selected_idx = stability_scores > 0.6
print(f'Selected {selected_idx.sum()} features (stability >0.6)')
```

**QC Checkpoint 2:**
- Selected features: 5-200 range
- Too few (<5): lower threshold, increase iterations
- Too many (>200): increase threshold, add pre-filtering

## Step 3: Leakage-Safe Performance Estimation

**Goal:** Estimate performance without the selection-before-CV leakage that inflates AUC toward 1.0 even on noise.

**Approach:** The Step 2 selection produced the discovery panel (fit on all training data) -- that is fine for the final panel, but it must NOT be the data the performance number is computed on. Estimate performance with scaling and selection wrapped in a `Pipeline` so they re-fit inside each fold; for raw RNA-seq, do per-sample normalization outside the fold and gene scaling/selection inside it.

```python
from sklearn.model_selection import StratifiedGroupKFold, cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression

# Selection lives INSIDE the pipeline -> re-fit per fold, no leakage. Use the unscaled X_train.
pipe = Pipeline([
    ('scaler', StandardScaler()),
    ('select', SelectKBest(f_classif, k=min(50, X_train.shape[1]))),
    ('clf', LogisticRegression(max_iter=5000, class_weight='balanced')),
])
# Group-aware outer CV: pass groups_train so no subject spans a fold boundary.
outer_cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)
cv_scores = cross_val_score(pipe, X_train, y_train, groups=groups_train, cv=outer_cv, scoring='roc_auc')
print(f'Leakage-safe CV AUC: {cv_scores.mean():.3f} +/- {cv_scores.std():.3f}')
# If hyperparameters are tuned, wrap a GridSearchCV (inner group CV) as the pipeline's estimator
# and report the OUTER cross_val_score -- flat CV that both tunes and reports is optimistic (Varma & Simon 2006).
```

**QC Checkpoint 3:**
- AUC reported with its fold spread, not a bare number (small-n CV is high-variance)
- Confirm selection is inside the pipeline and folds are group-aware; selection-before-CV inflates AUC toward 1.0 even on noise
- For imbalanced data report AUPRC/MCC, not accuracy; check the model predicts biology not batch (machine-learning/omics-classifiers)

## Step 4: Model Interpretation

**Goal:** Audit what the final model keys on, not select biomarkers.

**Approach:** Fit the final model on the discovery panel, then compute interventional SHAP against a background and aggregate over modules to catch shortcut/batch learning.

```python
import shap
import numpy as np
from sklearn.ensemble import RandomForestClassifier

# Fit the FINAL model on the discovery panel for interpretation and deployment.
sel = X_train.columns[selected_idx]
clf = RandomForestClassifier(n_estimators=300, random_state=42, n_jobs=-1).fit(X_train[sel], y_train)

# Interventional SHAP ('what the model uses') needs a background; set feature_perturbation
# explicitly because the 0.47+ 'auto' default flips the estimand on whether data= is given.
background = shap.utils.sample(X_train[sel], 100)
explainer = shap.TreeExplainer(clf, data=background, feature_perturbation='interventional')
shap_values = explainer(X_test[sel])
# RF returns one output per class in shap 0.47+ (n_samples, n_features, n_classes); keep the positive class.
if shap_values.values.ndim == 3:
    shap_values = shap_values[:, :, 1]
mean_shap = np.abs(shap_values.values).mean(axis=0)
```

**QC Checkpoint 4:**
- SHAP is an audit, not a selection method: use it to confirm the model is not keying on batch/housekeeping shortcuts (machine-learning/prediction-explanation)
- Aggregate SHAP over co-expression modules before ranking; within-module order is not a finding
- SHAP directions should be biologically plausible; treat top-SHAP genes as hypotheses, not a validated panel

## Step 5: Final Validation -- Discrimination AND Calibration

**Goal:** Report honest held-out performance, including calibration when risks will be used.

**Approach:** Report discrimination with an interval, but if the panel will produce risk estimates, also check calibration: AUC is invariant to any monotone transform of the score, so a high AUC says nothing about whether the probabilities are honest (machine-learning/model-validation). External validation on an independent cohort is the real bar.

```python
from sklearn.metrics import roc_auc_score
from sklearn.metrics import brier_score_loss
import numpy as np

y_prob = clf.predict_proba(X_test[sel])[:, 1]
test_auc = roc_auc_score(y_test, y_prob)

# Bootstrap CI for AUC (1000 resamples). Skip single-class resamples (roc_auc_score is nan there).
boot = []
for _ in range(1000):
    i = np.random.choice(len(y_test), len(y_test), replace=True)
    if len(np.unique(y_test[i])) == 2:
        boot.append(roc_auc_score(y_test[i], y_prob[i]))
ci_lower, ci_upper = np.percentile(boot, [2.5, 97.5])
print(f'Hold-out AUC: {test_auc:.3f}  95% CI [{ci_lower:.3f}, {ci_upper:.3f}]')
print(f'Brier score (calibration + refinement): {brier_score_loss(y_test, y_prob):.3f}')   # y must be 0/1-encoded; string labels raise unless pos_label is passed
# If risks will be used, recalibrate on a disjoint fold and report a reliability curve
# (machine-learning/model-validation); do not resample for imbalance -- it breaks calibration.
```

## Parameter Recommendations

| Step | Parameter | Recommendation |
|------|-----------|----------------|
| Split | n_splits (StratifiedGroupKFold) | 5 -> ~0.2 held out; lower n_splits for a larger test fraction |
| Boruta | max_iter | 100 (sufficient), 200 if tentative features |
| LASSO | n_bootstrap | 100 (quick), 500 for publication |
| LASSO | stability_threshold | 0.6 (standard), 0.8 for strict |
| Leakage-safe CV | folds | 5 (standard), 10 for small datasets; selection inside each fold |
| RF | n_estimators | 100-500 |
| XGBoost | learning_rate | 0.1 (conservative) |

## Common Errors

The leakage seams first (each silently inflates performance and a held-out set cannot detect it), then operational issues.

| Symptom | Cause | Fix |
|---------|-------|-----|
| Near-perfect CV AUC that collapses on external data | Features selected on the full dataset, then only the classifier CV'd | Wrap selection INSIDE the pipeline so it re-fits per fold (Ambroise & McLachlan 2002) |
| Optimistic AUC despite in-fold selection | Scaler/ComBat/PCA/imputation fit on all data before the split | Fit every data-dependent transform inside the fold (Pipeline) |
| Great CV, poor real-world performance | Repeated subjects (biopsies/longitudinal/replicates) split across train/test | Split by subject with StratifiedGroupKFold; the unit is the donor, not the sample |
| Reported AUC higher than any real fold | Same CV used to tune hyperparameters AND report | Nest: inner CV tunes, outer CV reports (Varma & Simon 2006) |
| Good AUC but risk estimates are miscalibrated | Resampling (SMOTE/undersampling) for imbalance, or AUC used as the only metric | Report AUPRC/MCC + Brier; recalibrate on a disjoint fold; do not resample-then-report calibration |
| No features selected | Too strict threshold | Lower stability threshold, increase iterations |
| Too many features (>200) | Noisy data | Add pre-filtering, increase regularization |
| Low CV AUC (<0.6) | No signal, low power | Check data quality, add samples |
| High variance across folds | Small sample size | Repeated stratified k-fold with an interval (LOOCV is degenerate for AUC) |
| SHAP features differ from selected | Correlated features split credit; attribution describes the model | Aggregate over modules; do not expect SHAP to match selection |

## Export Results

```python
import pandas as pd
import joblib

# Save biomarker panel
feature_names = X_train.columns[selected_idx].tolist()
pd.DataFrame({'feature': feature_names}).to_csv('biomarker_panel.csv', index=False)

# Save model and scaler for deployment
joblib.dump(clf, 'biomarker_classifier.joblib')
joblib.dump(scaler, 'feature_scaler.joblib')
```

## Related Skills

- database-access/geo-data - Public expression cohorts for validation sets
- database-access/sra-data - Pull raw FASTQ for re-quantified validation cohorts
- database-access/uniprot-access - Protein-level features (sequence, GO terms, PTMs) for protein biomarkers
- machine-learning/biomarker-discovery - Detailed feature selection methods
- machine-learning/model-validation - Nested CV implementation details
- machine-learning/omics-classifiers - Classifier options and tuning
- machine-learning/prediction-explanation - SHAP and LIME interpretation
- differential-expression/de-results - Pre-filter with DE genes
- pathway-analysis/go-enrichment - Functional enrichment of biomarkers

## References

- Ambroise C, McLachlan GJ (2002) Selection bias in gene extraction on the basis of microarray gene-expression data. *PNAS* 99:6562-6566. DOI 10.1073/pnas.102102699. (feature selection must be inside the CV fold.)
- Varma S, Simon R (2006) Bias in error estimation when using cross-validation for model selection. *BMC Bioinformatics* 7:91. DOI 10.1186/1471-2105-7-91. (nested CV for unbiased performance.)
- Whalen S, Schreiber J, Noble WS, Pollard KS (2022) Navigating the pitfalls of applying machine learning in genomics. *Nature Reviews Genetics* 23:169-181. DOI 10.1038/s41576-021-00434-9. (genomics-specific leakage and distribution-shift pitfalls.)
- Kapoor S, Narayanan A (2023) Leakage and the reproducibility crisis in machine-learning-based science. *Patterns* 4:100804. DOI 10.1016/j.patter.2023.100804. (a taxonomy of leakage, including group leakage.)
