---
name: bio-machine-learning-omics-classifiers
description: Builds classification models for omics data using RandomForest, XGBoost, and logistic regression with sklearn-compatible APIs. Includes proper preprocessing and evaluation metrics for biomarker classifiers. Use when building diagnostic or prognostic classifiers from expression or variant data.
tool_type: python
primary_tool: sklearn
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, pandas 2.2+, scikit-learn 1.4+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Classification Models for Omics Data

**"Build a classifier from my gene expression data"** → Train RandomForest, XGBoost, or logistic regression models on omics features with proper preprocessing and evaluation metrics.
- Python: `sklearn.ensemble.RandomForestClassifier()`, `xgboost.XGBClassifier()`

## Core Workflow

**Goal:** Train a classification model on omics data and evaluate its predictive performance.

**Approach:** Build a scaled pipeline with a Random Forest classifier, fit on training data, and assess with ROC-AUC on held-out test data.

```python
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, roc_auc_score, roc_curve
import matplotlib.pyplot as plt

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, stratify=y, random_state=42)

pipe = Pipeline([
    ('scaler', StandardScaler()),
    ('clf', RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1))
])
pipe.fit(X_train, y_train)

y_pred = pipe.predict(X_test)
y_prob = pipe.predict_proba(X_test)[:, 1]
print(classification_report(y_test, y_pred))
print(f'ROC-AUC: {roc_auc_score(y_test, y_prob):.3f}')
```

## XGBoost Classifier

**Goal:** Train a gradient-boosted tree classifier using the sklearn-compatible XGBoost API.

**Approach:** Configure XGBClassifier with proper parameter names (avoiding deprecated aliases) and wrap in a scaling pipeline.

```python
from xgboost import XGBClassifier

# Use sklearn-compatible API with proper parameters (avoid deprecated seed, nthread)
xgb = XGBClassifier(
    n_estimators=100,
    max_depth=6,
    learning_rate=0.1,
    random_state=42,  # NOT seed
    n_jobs=-1,        # NOT nthread
    eval_metric='logloss'
)

pipe = Pipeline([('scaler', StandardScaler()), ('clf', xgb)])
pipe.fit(X_train, y_train)
```

## Logistic Regression with Regularization

**Goal:** Build an interpretable linear classifier that simultaneously selects sparse biomarker features.

**Approach:** Use L1-regularized logistic regression with built-in cross-validation for penalty selection, then extract nonzero coefficients as selected features.

```python
from sklearn.linear_model import LogisticRegressionCV

# L1 for sparse biomarkers, L2 for correlated features, elasticnet for mixed
logit = LogisticRegressionCV(
    Cs=10,
    cv=5,
    penalty='l1',
    solver='saga',
    max_iter=1000,
    random_state=42
)
pipe = Pipeline([('scaler', StandardScaler()), ('clf', logit)])
pipe.fit(X_train, y_train)

# Get selected features (nonzero coefficients)
feature_mask = logit.coef_[0] != 0
selected = X.columns[feature_mask]
```

## ROC Curve Visualization

**Goal:** Generate a publication-quality ROC curve showing classifier discrimination ability.

**Approach:** Compute false/true positive rates from predicted probabilities and plot with AUC annotation.

```python
fpr, tpr, _ = roc_curve(y_test, y_prob)
auc = roc_auc_score(y_test, y_prob)

plt.figure(figsize=(6, 6))
plt.plot(fpr, tpr, label=f'ROC (AUC = {auc:.3f})')
plt.plot([0, 1], [0, 1], 'k--')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend()
plt.savefig('roc_curve.png', dpi=150)
```

## Multi-class Classification

**Goal:** Handle classification tasks with more than two classes while addressing class imbalance.

**Approach:** Encode labels numerically and use balanced class weights to upweight underrepresented classes during training.

```python
from sklearn.metrics import classification_report
from sklearn.preprocessing import LabelEncoder

le = LabelEncoder()
y_encoded = le.fit_transform(y)

# Use class_weight for imbalanced data
rf = RandomForestClassifier(n_estimators=100, class_weight='balanced', random_state=42)
```

## Feature Importance from Trees

**Goal:** Rank features by their contribution to tree-based classifier predictions.

**Approach:** Extract Gini importances from a fitted Random Forest and sort to identify top contributing features.

```python
import pandas as pd

importances = pipe.named_steps['clf'].feature_importances_
feature_imp = pd.DataFrame({'feature': X.columns, 'importance': importances})
feature_imp = feature_imp.sort_values('importance', ascending=False).head(20)
```

## Preprocessing Guidelines

| Data Type | Scaler | Notes |
|-----------|--------|-------|
| Log-counts (RNA-seq) | StandardScaler | Assumes ~normal after log |
| TPM/FPKM | StandardScaler | Gene-wise centering |
| Raw counts | None | Tree models handle counts |
| Mixed features | ColumnTransformer | Different scalers per type |

## Related Skills

- machine-learning/model-validation - Proper model evaluation
- machine-learning/prediction-explanation - Explain predictions with SHAP
- machine-learning/biomarker-discovery - Reduce features before modeling
