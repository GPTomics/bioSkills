'''Linear vs RF vs XGBoost on omics-like data: discrimination AND calibration.

Runs end-to-end on synthetic data with a purely LINEAR p>>n signal (no
interactions). Shows the two dossier points: a regularized linear model is
competitive with (here, beats) tree ensembles when the signal is linear and p>>n;
and tree ensembles are often miscalibrated -- a random forest bounds its votes
away from 0/1 while this XGBoost (log-loss, many rounds) pushes to the extremes --
so report Brier, not just AUC.
'''
# Reference: numpy 1.26+, scikit-learn 1.4+, xgboost 2.0+ | Verify API if version differs

import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegressionCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, brier_score_loss
from xgboost import XGBClassifier

rng = np.random.default_rng(0)
n, p, k = 300, 1500, 25                                # p >> n; 25 truly informative, linear
X = rng.normal(size=(n, p))
beta = np.zeros(p); beta[:k] = rng.normal(scale=1.2, size=k)
logit = X @ beta
y = (rng.random(n) < 1.0 / (1.0 + np.exp(-logit))).astype(int)          # linear generative model
X_tr, X_te, y_tr, y_te = train_test_split(X, y, test_size=0.4, stratify=y, random_state=0)

models = {
    # Genuine elastic-net (saga + l1_ratio); default-L2 leaves linear signal on the table.
    'elastic-net logistic': Pipeline([('s', StandardScaler()),
                                       ('c', LogisticRegressionCV(penalty='elasticnet', solver='saga',
                                                                  l1_ratios=[0.5], Cs=10, cv=5, max_iter=5000))]),
    'random forest': RandomForestClassifier(n_estimators=400, min_samples_leaf=3, random_state=0),
    'xgboost': XGBClassifier(n_estimators=300, learning_rate=0.05, max_depth=4,
                             subsample=0.8, colsample_bytree=0.5, eval_metric='logloss', random_state=0),
}

print(f'{"model":24s} {"AUC":>6s} {"Brier":>7s}')
for name, m in models.items():
    m.fit(X_tr, y_tr)
    pr = m.predict_proba(X_te)[:, 1]
    # AUC (ranking) can be similar across models while Brier (calibration) differs.
    print(f'{name:24s} {roc_auc_score(y_te, pr):6.3f} {brier_score_loss(y_te, pr):7.3f}')

# Tree probability spread vs logistic: RF compresses toward 0.5, boosting toward the extremes.
rf_p = models['random forest'].predict_proba(X_te)[:, 1]
xgb_p = models['xgboost'].predict_proba(X_te)[:, 1]
print(f'\nRF prob range:  [{rf_p.min():.2f}, {rf_p.max():.2f}]  (compressed toward 0.5)')
print(f'XGB prob range: [{xgb_p.min():.2f}, {xgb_p.max():.2f}]  (pushed toward extremes)')
