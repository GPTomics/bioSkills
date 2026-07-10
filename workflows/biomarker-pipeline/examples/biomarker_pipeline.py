'''End-to-end biomarker discovery pipeline'''
# Reference: matplotlib 3.8+, numpy 1.26+, pandas 2.2+, scikit-learn 1.4+, boruta 0.4+, shap 0.47+, joblib 1.3+ | Verify API if version differs

import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedGroupKFold, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, classification_report, brier_score_loss
from sklearn.calibration import calibration_curve
from boruta import BorutaPy
import shap
import matplotlib.pyplot as plt

# Load data
# Example data: Use GEO datasets (e.g., GSE37418) or Bioconductor's curatedOvarianData
expr = pd.read_csv('expression.csv', index_col=0)
meta = pd.read_csv('metadata.csv', index_col=0)
X = expr.T  # transpose to samples x genes
# y must be 0/1: brier_score_loss and calibration_curve raise on string labels. Encode the DISEASE
# class as 1 explicitly -- LabelEncoder assigns 1 alphabetically, which would make 'control' positive
# for a case/control column and silently flip what y_prob, the SHAP slice, and Brier describe.
POSITIVE_CLASS = 'disease'
y = (meta.loc[X.index, 'condition'].values == POSITIVE_CLASS).astype(int)
# The independent unit is the SUBJECT, not the sample. If metadata lacks a subject id assume
# one sample per subject (groups = X.index); otherwise this MUST be the subject key or repeated
# subjects leak across the split and a held-out set cannot detect it.
groups = meta['subject_id'].loc[X.index].values if 'subject_id' in meta.columns else np.asarray(X.index)

print(f'Data: {X.shape[0]} samples, {X.shape[1]} features')
print(f'Classes: {np.unique(y, return_counts=True)}')

# Step 1: Group- AND class-aware hold-out -- one StratifiedGroupKFold fold (~1/5) as the test set,
# so no subject spans train/test. train_test_split(stratify=y) alone would leak repeated subjects.
sgkf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)
train_idx, test_idx = next(sgkf.split(X, y, groups))
X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
y_train, y_test = y[train_idx], y[test_idx]
groups_train = groups[train_idx]
assert not (set(groups[train_idx]) & set(groups[test_idx])), 'subject leaked across train/test'
print(f'Train: {len(y_train)}, Test: {len(y_test)}')

# Step 2: Discovery panel -- fit scaler + Boruta on ALL training data.
# This panel is the deliverable; the honest performance estimate comes from Step 3, not here.
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# max_depth=5: Shallow trees for stable importances across Boruta iterations
rf_selector = RandomForestClassifier(n_estimators=100, max_depth=5, n_jobs=-1, random_state=42)
# max_iter=100: Usually sufficient; increase to 200 if many tentative features remain
# n_estimators='auto': BorutaPy sets tree count from int((sqrt(2*n_features)/max_depth)*100)
boruta = BorutaPy(rf_selector, n_estimators='auto', max_iter=100, random_state=42, verbose=0)
boruta.fit(X_train_scaled, y_train)

selected_features = X_train.columns[boruta.support_].tolist()
print(f'Selected {len(selected_features)} features')

# QC: Check feature count is in reasonable range (5-200)
if len(selected_features) < 5:
    print('WARNING: Few features selected. Consider lowering threshold or increasing max_iter.')
elif len(selected_features) > 200:
    print('WARNING: Many features selected. Consider stricter pre-filtering.')

X_train_sel = X_train_scaled[:, boruta.support_]
X_test_sel = X_test_scaled[:, boruta.support_]

# Step 3: Leakage-safe performance estimate -- selection runs INSIDE each fold.
# Re-selecting per fold is the only honest estimate; cross_val_score on the Step-2
# panel would leak the whole training set into selection and inflate AUC toward 1.0.
outer_cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)   # group-aware folds
leakage_safe = Pipeline([
    ('scale', StandardScaler()),
    ('select', BorutaPy(RandomForestClassifier(n_estimators=100, max_depth=5, n_jobs=-1, random_state=42),
                        n_estimators='auto', max_iter=100, random_state=42, verbose=0)),
    ('clf', RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1)),
])
# error_score=np.nan so a fold where Boruta selects 0 features (common on weak-signal omics)
# yields nan for that fold instead of aborting the whole CV.
cv_scores = cross_val_score(leakage_safe, X_train.values, y_train, groups=groups_train,
                            cv=outer_cv, scoring='roc_auc', error_score=np.nan)

print(f'Leakage-safe CV AUC: {np.nanmean(cv_scores):.3f} +/- {np.nanstd(cv_scores):.3f}')

# QC: Check AUC and variance
if np.nanmean(cv_scores) < 0.7:
    print('WARNING: Low AUC. Check data quality or add samples.')
if np.nanstd(cv_scores) > 0.1:
    # Repeated stratified group k-fold for a tighter estimate -- NOT LOOCV (single-sample folds
    # make per-fold roc_auc undefined; see SKILL.md).
    print('WARNING: High fold variance. Use repeated stratified group k-fold and report an interval.')

# Step 4: Refit the panel classifier and audit it with interventional SHAP on HELD-OUT data.
clf = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1)
clf.fit(X_train_sel, y_train)

# Set feature_perturbation explicitly: 'auto' silently flips estimand by shap version.
# Interventional (marginal) attributions need a background sample; audit on the test fold.
background = shap.utils.sample(X_train_sel, 100, random_state=42)
explainer = shap.TreeExplainer(clf, data=background, feature_perturbation='interventional')
shap_values = explainer(X_test_sel)

# Binary classifier returns one output per class; keep the positive class for plotting.
if shap_values.values.ndim == 3:
    shap_values = shap_values[:, :, 1]

# Beeswarm plot: shows importance AND direction
# max_display=20: Top 20 features for readability
shap.plots.beeswarm(shap_values, max_display=20, show=False)
plt.tight_layout()
plt.savefig('shap_beeswarm.png', dpi=150, bbox_inches='tight')
plt.close()
print('Saved SHAP beeswarm plot')

# Extract top SHAP features for QC comparison
mean_shap = np.abs(shap_values.values).mean(axis=0)
top_shap_idx = np.argsort(mean_shap)[-20:]
shap_feature_df = pd.DataFrame({
    'feature': [selected_features[i] for i in top_shap_idx],
    'mean_shap': mean_shap[top_shap_idx]
}).sort_values('mean_shap', ascending=False)
shap_feature_df.to_csv('shap_top_features.csv', index=False)

# Step 5: Validate on hold-out test set
y_prob = clf.predict_proba(X_test_sel)[:, 1]
test_auc = roc_auc_score(y_test, y_prob)
print(f'Hold-out test AUC: {test_auc:.3f}')

# Bootstrap CI for AUC (n_bootstrap=1000 for publication-quality intervals).
# Skip single-class resamples: roc_auc_score returns nan on a resample with one class (likely on a
# small/imbalanced hold-out), which would poison the percentile CI.
n_bootstrap = 1000
boot_aucs = []
for i in range(n_bootstrap):
    idx = np.random.choice(len(y_test), size=len(y_test), replace=True)
    if len(np.unique(y_test[idx])) < 2:
        continue
    boot_aucs.append(roc_auc_score(y_test[idx], y_prob[idx]))
ci_lower, ci_upper = np.percentile(boot_aucs, [2.5, 97.5])   # 95% CI
print(f'95% CI: [{ci_lower:.3f}, {ci_upper:.3f}]')

# Discrimination is not enough -- report CALIBRATION too (governing principle #4).
brier = brier_score_loss(y_test, y_prob)
print(f'Brier score (lower is better): {brier:.3f}')
frac_pos, mean_pred = calibration_curve(y_test, y_prob, n_bins=5, strategy='quantile')
print('Reliability (mean predicted vs observed frequency):')
for mp, fp in zip(mean_pred, frac_pos):
    print(f'  pred {mp:.2f} -> observed {fp:.2f}')

# Classification report
print('\nClassification Report:')
print(classification_report(y_test, clf.predict(X_test_sel), target_names=[f'not_{POSITIVE_CLASS}', POSITIVE_CLASS]))

# Export results
pd.DataFrame({'feature': selected_features}).to_csv('biomarker_panel.csv', index=False)
print(f'\nExported {len(selected_features)} biomarkers to biomarker_panel.csv')

# Optional: Save model for deployment
import joblib
joblib.dump(clf, 'biomarker_classifier.joblib')
joblib.dump(scaler, 'feature_scaler.joblib')
print('Saved classifier and scaler')
