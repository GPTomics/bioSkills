'''The SHAP credit split between correlated genes is mode-dependent (not biology).

Runs end-to-end on synthetic data. Two co-expressed genes both carry the signal;
TreeSHAP splits credit between them differently under tree_path_dependent
(conditional) vs interventional (marginal) conditioning. The module TOTAL is
stable, but the within-module order is an artifact of the chosen estimand --
which is why you aggregate over modules before ranking.
'''
# Reference: numpy 1.26+, scikit-learn 1.4+, shap 0.44+ | Verify API if version differs

import numpy as np
import shap
from sklearn.ensemble import RandomForestClassifier

rng = np.random.default_rng(0)
n = 800
A = rng.normal(size=n)
B = A + rng.normal(scale=0.8, size=n)               # B co-expressed with A but noisier (corr ~0.78)
noise = rng.normal(size=(n, 8))
X = np.column_stack([A, B, noise])
y = (A + rng.normal(scale=0.3, size=n) > 0).astype(int)    # label depends on A; the model still uses B as a noisy proxy

model = RandomForestClassifier(n_estimators=300, random_state=0).fit(X, y)
print(f'corr(geneA, geneB) = {np.corrcoef(A, B)[0, 1]:.2f}')

def module_split(mode, data=None):
    expl = shap.TreeExplainer(model, data=data, feature_perturbation=mode)
    vals = expl.shap_values(X[:200])
    vals = vals[1] if isinstance(vals, list) else (vals[..., 1] if vals.ndim == 3 else vals)
    m = np.abs(vals).mean(axis=0)
    return m[0], m[1], m[0] + m[1]

pa, pb, ptot = module_split('tree_path_dependent')
ia, ib, itot = module_split('interventional', data=shap.utils.sample(X, 100))
print(f'{"mode":20s} {"geneA":>7s} {"geneB":>7s} {"A+B":>7s}')
print(f'{"path_dependent":20s} {pa:7.3f} {pb:7.3f} {ptot:7.3f}')
print(f'{"interventional":20s} {ia:7.3f} {ib:7.3f} {itot:7.3f}')
print('\nThe A-vs-B split differs by mode; the module total is the stable quantity.')
print('Lesson: aggregate |SHAP| over co-expression modules before ranking genes.')
