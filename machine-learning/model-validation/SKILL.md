---
name: bio-machine-learning-model-validation
description: Implements nested cross-validation and stratified splits for unbiased model evaluation on biomedical datasets. Prevents data leakage and overfitting in biomarker discovery. Use when validating classifiers or optimizing hyperparameters on omics data.
tool_type: python
primary_tool: sklearn
---

## Version Compatibility

Reference examples tested with: numpy 1.26+, scikit-learn 1.4+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Cross-Validation for Biomedical Data

**"Properly validate my omics classifier"** → Use nested cross-validation with stratified splits to get unbiased performance estimates while tuning hyperparameters on small biomedical datasets.
- Python: `sklearn.model_selection.cross_val_score()` with `StratifiedKFold` inner/outer loops

## Why Nested CV Matters

Simple train/test splits overestimate performance on small omics datasets. Nested CV provides unbiased estimates by separating hyperparameter tuning from performance evaluation.

## Nested Cross-Validation

**Goal:** Obtain unbiased performance estimates by separating hyperparameter tuning from evaluation.

**Approach:** Use an outer CV loop for scoring and an inner CV loop for grid search, preventing information leakage between tuning and evaluation.

```python
from sklearn.model_selection import cross_val_score, StratifiedKFold, GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
import numpy as np

pipe = Pipeline([
    ('scaler', StandardScaler()),
    ('clf', RandomForestClassifier(random_state=42))
])

param_grid = {
    'clf__n_estimators': [50, 100, 200],
    'clf__max_depth': [5, 10, None]
}

# Outer CV: performance estimation (5 folds)
# Inner CV: hyperparameter tuning (3 folds)
outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
inner_cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=42)

nested_scores = []
for train_idx, test_idx in outer_cv.split(X, y):
    X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]

    grid = GridSearchCV(pipe, param_grid, cv=inner_cv, scoring='roc_auc', n_jobs=-1)
    grid.fit(X_train, y_train)
    score = grid.score(X_test, y_test)
    nested_scores.append(score)

print(f'Nested CV AUC: {np.mean(nested_scores):.3f} +/- {np.std(nested_scores):.3f}')
```

## Stratified K-Fold

**Goal:** Evaluate model performance while preserving class proportions in each fold.

**Approach:** Split data into stratified folds and compute cross-validated scores to account for class imbalance.

```python
from sklearn.model_selection import StratifiedKFold, cross_val_score

# Always stratify for class imbalance
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
scores = cross_val_score(pipe, X, y, cv=cv, scoring='roc_auc')
print(f'CV AUC: {scores.mean():.3f} +/- {scores.std():.3f}')
```

## Repeated Stratified K-Fold

**Goal:** Produce more stable performance estimates by averaging across multiple CV repetitions.

**Approach:** Repeat stratified K-fold splitting with different random seeds and aggregate scores across all iterations.

```python
from sklearn.model_selection import RepeatedStratifiedKFold

# More robust estimate with multiple repeats
cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=42)
scores = cross_val_score(pipe, X, y, cv=cv, scoring='roc_auc')
print(f'Repeated CV AUC: {scores.mean():.3f} +/- {scores.std():.3f}')
```

## Leave-One-Out (Small Datasets)

**Goal:** Maximize training data when sample size is very small (n < 30).

**Approach:** Hold out one sample at a time for testing and train on all remaining samples, then aggregate predictions.

```python
from sklearn.model_selection import LeaveOneOut, cross_val_predict

# Use for very small datasets (n < 30)
loo = LeaveOneOut()
y_pred = cross_val_predict(pipe, X, y, cv=loo, method='predict_proba')[:, 1]
auc = roc_auc_score(y, y_pred)
print(f'LOO AUC: {auc:.3f}')
```

## Group-Aware Splits

**Goal:** Prevent data leakage when samples from the same patient or batch are correlated.

**Approach:** Use group-aware splitting to ensure all samples from a single group stay in the same fold.

```python
from sklearn.model_selection import GroupKFold, LeaveOneGroupOut

# When samples from same patient/batch must stay together
groups = meta['patient_id'].values
group_cv = GroupKFold(n_splits=5)
scores = cross_val_score(pipe, X, y, cv=group_cv, groups=groups, scoring='roc_auc')
```

## CV Strategy Selection

| Dataset Size | Strategy | Notes |
|--------------|----------|-------|
| n > 100 | StratifiedKFold(5) | Standard choice |
| n = 50-100 | StratifiedKFold(10) | More train data per fold |
| n < 30 | LeaveOneOut | Maximum train data |
| Repeated measures | GroupKFold | Keep patients together |
| High variance | RepeatedStratifiedKFold | More stable estimates |

## Avoiding Data Leakage

**Goal:** Ensure feature selection does not use test-fold information, which inflates performance estimates.

**Approach:** Embed feature selection inside a pipeline so it executes independently within each CV fold.

```python
# WRONG: Feature selection before CV
# selected = SelectKBest(k=100).fit_transform(X, y)  # Leaks info!
# scores = cross_val_score(clf, selected, y, cv=cv)

# CORRECT: Feature selection inside CV
from sklearn.feature_selection import SelectKBest

pipe = Pipeline([
    ('scaler', StandardScaler()),
    ('select', SelectKBest(k=100)),  # Done per fold
    ('clf', RandomForestClassifier())
])
scores = cross_val_score(pipe, X, y, cv=cv, scoring='roc_auc')
```

## Related Skills

- machine-learning/omics-classifiers - Model training
- experimental-design/multiple-testing - Multiple hypothesis correction
- machine-learning/biomarker-discovery - Feature selection within CV
