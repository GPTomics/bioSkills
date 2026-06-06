'''Nested CV with discrimination AND calibration, on synthetic omics-like data.

Demonstrates the two things this skill insists on, end-to-end on synthetic data:
1. Nested CV (selection + scaling inside the fold) gives an honest AUC, while
   selecting features on the full matrix first inflates it.
2. AUC and calibration are orthogonal -- a model can discriminate yet be
   miscalibrated; report a reliability curve and the Brier score, then recalibrate.
'''
# Reference: numpy 1.26+, scikit-learn 1.4+ | Verify API if version differs

import numpy as np
from sklearn.datasets import make_classification
from sklearn.model_selection import StratifiedKFold, GridSearchCV, cross_val_score, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.calibration import calibration_curve, CalibratedClassifierCV
from sklearn.metrics import roc_auc_score, brier_score_loss

X, y = make_classification(n_samples=300, n_features=2000, n_informative=15,
                           weights=[0.85, 0.15], random_state=0)        # p>>n, imbalanced

# --- 1. Nested CV vs leaky selection ---
pipe = Pipeline([('scaler', StandardScaler()),
                 ('select', SelectKBest(f_classif)),
                 ('clf', LogisticRegression(max_iter=5000))])
grid = {'select__k': [20, 100], 'clf__C': [0.01, 0.1, 1]}
inner = StratifiedKFold(5, shuffle=True, random_state=0)
outer = StratifiedKFold(5, shuffle=True, random_state=1)
nested = cross_val_score(GridSearchCV(pipe, grid, cv=inner, scoring='roc_auc'),
                         X, y, cv=outer, scoring='roc_auc')

leaky_k = SelectKBest(f_classif, k=100).fit(X, y).get_support(indices=True)   # selected on ALL data
leaky = cross_val_score(LogisticRegression(max_iter=5000), X[:, leaky_k], y, cv=outer, scoring='roc_auc')
print(f'Nested-safe AUC: {nested.mean():.3f} +/- {nested.std():.3f}')
print(f'Leaky (selection-before-CV) AUC: {leaky.mean():.3f}   <- optimistic')

# --- 2. Discrimination vs calibration: oversampling the minority in training
# --- inflates predicted risk (the imbalance-correction trap) for no AUC gain. ---
Xc, yc = make_classification(n_samples=2000, n_features=20, n_informative=12,
                             weights=[0.9, 0.1], random_state=1)
X_tr, X_te, y_tr, y_te = train_test_split(Xc, yc, test_size=0.4, stratify=yc, random_state=0)
X_tr, X_cal, y_tr, y_cal = train_test_split(X_tr, y_tr, test_size=0.33, stratify=y_tr, random_state=0)

pos = np.where(y_tr == 1)[0]
balanced = np.concatenate([np.arange(len(y_tr)), np.repeat(pos, 8)])    # crude oversample to ~50/50
clf = LogisticRegression(max_iter=5000).fit(X_tr[balanced], y_tr[balanced])
p_raw = clf.predict_proba(X_te)[:, 1]
prevalence = y_te.mean()
print(f'\nLogReg on oversampled data AUC: {roc_auc_score(y_te, p_raw):.3f}   Brier (raw): {brier_score_loss(y_te, p_raw):.3f}')
print(f'Mean predicted risk {p_raw.mean():.2f} vs true prevalence {prevalence:.2f}  <- risk inflated by oversampling')

# Recalibrate on a disjoint fold (cv='prefit' deprecated >=1.6; FrozenEstimator on newer sklearn).
try:
    from sklearn.frozen import FrozenEstimator
    cal = CalibratedClassifierCV(FrozenEstimator(clf), method='isotonic').fit(X_cal, y_cal)
except ImportError:
    cal = CalibratedClassifierCV(clf, method='isotonic', cv='prefit').fit(X_cal, y_cal)
p_cal = cal.predict_proba(X_te)[:, 1]
print(f'Calibrated Brier: {brier_score_loss(y_te, p_cal):.3f}   (AUC ~unchanged: {roc_auc_score(y_te, p_cal):.3f}; isotonic can nudge it via ties)')

prob_true, prob_pred = calibration_curve(y_te, p_cal, n_bins=5, strategy='quantile')
print('Reliability (pred -> observed):', list(zip(prob_pred.round(2), prob_true.round(2))))
