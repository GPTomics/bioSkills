'''L1-regularized logistic regression for sparse biomarker selection'''

import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import classification_report, roc_auc_score

expr = pd.read_csv('expression.csv', index_col=0)
meta = pd.read_csv('metadata.csv', index_col=0)

X = expr.T
y = meta.loc[X.index, 'condition'].values

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, stratify=y, random_state=42)

# L1 (lasso): Sparse solutions, good for biomarker discovery
# L2 (ridge): Better for correlated features, keeps all
# elasticnet: Mix of both (l1_ratio controls balance)
logit = LogisticRegressionCV(
    Cs=10,         # 10 regularization strengths to try
    cv=5,          # 5-fold CV for tuning
    penalty='l1',  # L1 for sparsity
    solver='saga', # Required for L1
    max_iter=1000, # Increase for convergence
    random_state=42
)

pipe = Pipeline([('scaler', StandardScaler()), ('clf', logit)])
pipe.fit(X_train, y_train)

y_pred = pipe.predict(X_test)
y_prob = pipe.predict_proba(X_test)[:, 1]
print(classification_report(y_test, y_pred))
print(f'ROC-AUC: {roc_auc_score(y_test, y_prob):.3f}')
print(f'Best C: {logit.C_[0]:.4f}')

coefs = pd.Series(logit.coef_[0], index=X.columns)
nonzero = coefs[coefs != 0].sort_values(key=abs, ascending=False)
print(f'\nSelected {len(nonzero)} features with non-zero coefficients:')
nonzero.to_csv('logistic_biomarkers.csv')
nonzero.head(20)
