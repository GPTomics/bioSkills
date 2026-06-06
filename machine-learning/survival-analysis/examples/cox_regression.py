'''Predictive survival modeling: penalized Cox vs RSF, evaluated beyond the C-index.

Runs end-to-end on synthetic data. Fits an elastic-net Cox and a random survival
forest, then evaluates with Uno's IPCW C, time-dependent AUC, and integrated Brier
versus a Kaplan-Meier baseline -- because the C-index alone is censoring-dependent,
blind to calibration, and insensitive.
'''
# Reference: scikit-survival 0.22+, numpy 1.26+ | Verify API if version differs

import numpy as np
from sksurv.util import Surv
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import concordance_index_ipcw, cumulative_dynamic_auc, integrated_brier_score
from sklearn.model_selection import train_test_split

rng = np.random.default_rng(0)
n, p = 600, 40
X = rng.normal(size=(n, p))
beta = np.zeros(p); beta[:5] = 0.8                       # 5 prognostic features
risk = X @ beta
true_time = rng.exponential(np.exp(-risk))               # higher risk -> shorter time
censor = rng.exponential(1.5, n)
time = np.minimum(true_time, censor).astype(float)
event = (true_time <= censor)                            # bool event field is REQUIRED by sksurv

Xtr, Xte, etr, ete, ttr, tte = train_test_split(X, event, time, test_size=0.4, random_state=0)
y_tr = Surv.from_arrays(event=etr, time=ttr)             # structured array: (bool event, float time)
y_te = Surv.from_arrays(event=ete, time=tte)

cox = CoxnetSurvivalAnalysis(l1_ratio=0.9, alpha_min_ratio=0.01, fit_baseline_model=True).fit(Xtr, y_tr)
rsf = RandomSurvivalForest(n_estimators=300, min_samples_leaf=15, random_state=0).fit(Xtr, y_tr)

tau = np.quantile(tte[ete], 0.8)                         # truncate IPCW at a horizon inside follow-up
times = np.linspace(np.quantile(tte[ete], 0.1), tau, 12)

print(f'{"model":18s} {"UnoC":>6s} {"meanAUC":>8s} {"IBS":>6s}')
for name, m in [('elastic-net Cox', cox), ('random forest', rsf)]:
    risk_score = m.predict(Xte)                          # higher = higher risk, NOT a probability
    c = concordance_index_ipcw(y_tr, y_te, risk_score, tau=tau)[0]
    _, mean_auc = cumulative_dynamic_auc(y_tr, y_te, risk_score, times)
    surv = np.vstack([[fn(t) for t in times] for fn in m.predict_survival_function(Xte)])
    ibs = integrated_brier_score(y_tr, y_te, surv, times)
    print(f'{name:18s} {c:6.3f} {mean_auc:8.3f} {ibs:6.3f}')

# KM-only baseline IBS: a model whose IBS does not beat this has no predictive value.
km_surv = np.tile([np.mean(ttr[etr] > t) for t in times], (len(tte), 1))
print(f'\nKaplan-Meier baseline IBS (no covariates): {integrated_brier_score(y_tr, y_te, km_surv, times):.3f}')
