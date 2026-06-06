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

Reference examples tested with: numpy 1.26+, pandas 2.2+, scikit-learn 1.4+, shap 0.44+ (xgboost 2.0+ optional).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

scikit-learn drift: `CalibratedClassifierCV(cv='prefit')` deprecated in 1.6 (use `FrozenEstimator`); `LogisticRegression(penalty=)` deprecated in 1.8. XGBoost moved `early_stopping_rounds` to the constructor in 2.x. If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Biomarker Discovery Pipeline

**"Build a validated biomarker panel from my omics data"** -> Orchestrate feature selection, leakage-safe cross-validation, calibration, and SHAP interpretation to produce a robust, honestly-validated biomarker signature.

The one mistake that invalidates the whole pipeline: estimating performance after selecting features on the full training set. Selection is the dominant overfitting capacity in p>>n and gives near-perfect apparent accuracy on pure noise (machine-learning/model-validation). The discovery panel may be selected on all training data, but every reported performance number must come from a pipeline that re-runs selection inside each CV fold.

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

**Goal:** Load the matrix and hold out a test set before anything is fit.

**Approach:** Stratified split, then fit the scaler on training only; per-fold scaling is re-applied inside the CV pipeline in Step 3.

```python
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

expr = pd.read_csv('expression.csv', index_col=0)
meta = pd.read_csv('metadata.csv', index_col=0)

X = expr.T  # samples x genes
y = meta.loc[X.index, 'condition'].values

# test_size=0.2: Standard 80/20 split; use 0.3 for <100 samples
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, stratify=y, random_state=42
)

# Fit scaler on training only to prevent data leakage
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)
```

**QC Checkpoint 1:** Check class balance, sample counts per group
- Minimum 10 samples per class recommended
- Classes should be reasonably balanced (ratio <3:1)

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
from sklearn.model_selection import StratifiedKFold, cross_val_score
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
outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
cv_scores = cross_val_score(pipe, X_train, y_train, cv=outer_cv, scoring='roc_auc')
print(f'Leakage-safe CV AUC: {cv_scores.mean():.3f} +/- {cv_scores.std():.3f}')
```

**QC Checkpoint 3:**
- AUC reported with its fold spread, not a bare number (small-n CV is high-variance)
- Confirm selection is inside the pipeline; selection-before-CV inflates AUC toward 1.0 even on noise
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

# Bootstrap CI for AUC (1000 resamples for a publication-quality interval).
boot = [roc_auc_score(y_test[i], y_prob[i]) for i in
        (np.random.choice(len(y_test), len(y_test), replace=True) for _ in range(1000))]
ci_lower, ci_upper = np.percentile(boot, [2.5, 97.5])
print(f'Hold-out AUC: {test_auc:.3f}  95% CI [{ci_lower:.3f}, {ci_upper:.3f}]')
print(f'Brier score (calibration + refinement): {brier_score_loss(y_test, y_prob):.3f}')
# If risks will be used, recalibrate on a disjoint fold and report a reliability curve
# (machine-learning/model-validation); do not resample for imbalance -- it breaks calibration.
```

## Parameter Recommendations

| Step | Parameter | Recommendation |
|------|-----------|----------------|
| Split | test_size | 0.2 (standard), 0.3 for small datasets |
| Boruta | max_iter | 100 (sufficient), 200 if tentative features |
| LASSO | n_bootstrap | 100 (quick), 500 for publication |
| LASSO | stability_threshold | 0.6 (standard), 0.8 for strict |
| Leakage-safe CV | folds | 5 (standard), 10 for small datasets; selection inside each fold |
| RF | n_estimators | 100-500 |
| XGBoost | learning_rate | 0.1 (conservative) |

## Troubleshooting

| Issue | Likely Cause | Solution |
|-------|--------------|----------|
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
