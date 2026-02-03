'''Nested cross-validation for unbiased biomarker classifier evaluation'''

import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold, GridSearchCV, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, classification_report

expr = pd.read_csv('expression.csv', index_col=0)
meta = pd.read_csv('metadata.csv', index_col=0)

X = expr.T
y = meta.loc[X.index, 'condition'].values

pipe = Pipeline([
    ('scaler', StandardScaler()),
    ('clf', RandomForestClassifier(random_state=42, n_jobs=-1))
])

param_grid = {
    'clf__n_estimators': [50, 100, 200],
    'clf__max_depth': [5, 10, None],
    'clf__min_samples_leaf': [1, 3, 5]
}

# n_splits=5 outer: Standard for n>50; use 10 for smaller datasets
# n_splits=3 inner: Fewer folds for faster tuning, sufficient for hyperparameter selection
outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
inner_cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=42)

nested_scores = []
best_params = []

for fold, (train_idx, test_idx) in enumerate(outer_cv.split(X, y), 1):
    X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]

    grid = GridSearchCV(pipe, param_grid, cv=inner_cv, scoring='roc_auc', n_jobs=-1)
    grid.fit(X_train, y_train)

    y_prob = grid.predict_proba(X_test)[:, 1]
    score = roc_auc_score(y_test, y_prob)
    nested_scores.append(score)
    best_params.append(grid.best_params_)
    print(f'Fold {fold}: AUC = {score:.3f}, best params = {grid.best_params_}')

print(f'\nNested CV AUC: {np.mean(nested_scores):.3f} +/- {np.std(nested_scores):.3f}')

# Compare to non-nested (optimistically biased)
grid = GridSearchCV(pipe, param_grid, cv=5, scoring='roc_auc', n_jobs=-1)
grid.fit(X, y)
# cv_results_ best_score_ is optimistic due to selection bias
print(f'Non-nested CV AUC: {grid.best_score_:.3f} (biased estimate)')
print(f'Difference shows optimism: {grid.best_score_ - np.mean(nested_scores):.3f}')
