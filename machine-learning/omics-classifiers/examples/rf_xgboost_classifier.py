'''Train RandomForest and XGBoost classifiers on expression data'''

import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, roc_auc_score, roc_curve
from xgboost import XGBClassifier
import matplotlib.pyplot as plt

expr = pd.read_csv('expression.csv', index_col=0)  # genes x samples
meta = pd.read_csv('metadata.csv', index_col=0)

X = expr.T
y = meta.loc[X.index, 'condition'].values

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, stratify=y, random_state=42)

rf_pipe = Pipeline([
    ('scaler', StandardScaler()),
    ('clf', RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1))
])
rf_pipe.fit(X_train, y_train)

# n_estimators=100: XGBoost default; sufficient for most omics classification
# max_depth=6: XGBoost default; prevents overfitting when features >> samples
# learning_rate=0.1: Lower than default (0.3); slower learning reduces overfitting risk
# random_state=42: Reproducibility
xgb_pipe = Pipeline([
    ('scaler', StandardScaler()),
    ('clf', XGBClassifier(n_estimators=100, max_depth=6, learning_rate=0.1, random_state=42, n_jobs=-1, eval_metric='logloss'))
])
xgb_pipe.fit(X_train, y_train)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for ax, (name, pipe) in zip(axes, [('RandomForest', rf_pipe), ('XGBoost', xgb_pipe)]):
    y_prob = pipe.predict_proba(X_test)[:, 1]
    fpr, tpr, _ = roc_curve(y_test, y_prob)
    auc = roc_auc_score(y_test, y_prob)
    ax.plot(fpr, tpr, label=f'AUC = {auc:.3f}')
    ax.plot([0, 1], [0, 1], 'k--')
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title(name)
    ax.legend()
    print(f'\n{name}:')
    print(classification_report(y_test, pipe.predict(X_test)))

plt.tight_layout()
plt.savefig('roc_comparison.png', dpi=150)

importances = rf_pipe.named_steps['clf'].feature_importances_
top_features = pd.Series(importances, index=X.columns).nlargest(20)
top_features.to_csv('top_rf_features.csv')
print(f'\nTop 20 RF features saved')
