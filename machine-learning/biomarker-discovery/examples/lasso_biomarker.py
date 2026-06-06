'''Leakage-safe selection and stability, demonstrated on synthetic p>>n data.

Two teaching points run end-to-end on synthetic data:
1. Selecting features on the FULL matrix before CV reports high AUC even when the
   labels are pure noise; putting selection inside the CV fold does not.
2. A signature is reported with a stability index, not on its own.
'''
# Reference: numpy 1.26+, pandas 2.2+, scikit-learn 1.4+ | Verify API if version differs

import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.model_selection import cross_val_score, StratifiedKFold

rng = np.random.default_rng(0)
n, p = 80, 5000                                  # p >> n: 80 samples, 5000 'genes'
X = pd.DataFrame(rng.normal(size=(n, p)), columns=[f'g{i}' for i in range(p)])
y_noise = rng.integers(0, 2, n)                  # pure noise: no gene is truly associated

# WRONG: select top-20 on all data, then CV the classifier on those 20.
top20 = SelectKBest(f_classif, k=20).fit(X, y_noise).get_support(indices=True)
leaky = cross_val_score(LogisticRegression(max_iter=5000), X.iloc[:, top20], y_noise,
                        cv=StratifiedKFold(10, shuffle=True, random_state=0), scoring='roc_auc')

# RIGHT: selection inside the Pipeline, re-fit per fold.
pipe = Pipeline([('select', SelectKBest(f_classif, k=20)),
                 ('clf', LogisticRegression(max_iter=5000))])
honest = cross_val_score(pipe, X, y_noise, cv=StratifiedKFold(10, shuffle=True, random_state=0), scoring='roc_auc')

print(f'AUC on PURE NOISE, selection-before-CV (WRONG): {leaky.mean():.2f}')   # inflated, ~0.7+
print(f'AUC on PURE NOISE, selection-inside-CV (RIGHT): {honest.mean():.2f}')  # ~0.5 as it should be

# --- Stability on a separate matrix with real signal in the first 5 genes.
# --- Moderate p after a notional pre-filter so L1 can recover the signal. ---
ns, ps = 120, 300
Xs = pd.DataFrame(rng.normal(size=(ns, ps)), columns=[f'g{i}' for i in range(ps)])
beta = np.zeros(ps); beta[:5] = 2.5
logit = Xs.values @ beta + rng.normal(scale=1.0, size=ns)
ys = (logit > np.median(logit)).astype(int)

counts = np.zeros(ps); subsets = []
for _ in range(100):
    idx = rng.choice(ns, size=ns // 2, replace=False)          # n/2 subsampling
    mask = LogisticRegression(penalty='l1', solver='liblinear', C=1.0,
                              max_iter=2000).fit(Xs.iloc[idx], ys[idx]).coef_[0] != 0
    counts += mask; subsets.append(mask.astype(int))

Z = np.array(subsets); k = Z.sum(axis=1)
# Nogueira 2018 stability index: chance-corrected, 1 = identical, ~0 = random.
stability = 1 - Z.var(axis=0, ddof=1).mean() / ((k.mean() / ps) * (1 - k.mean() / ps))
stable = Xs.columns[counts / 100 > 0.6]                         # pi_thr = 0.6
print(f'\nStable features (>60%, true signal = g0..g4): {list(stable)}')
print(f'Nogueira stability index: {stability:.2f}')
