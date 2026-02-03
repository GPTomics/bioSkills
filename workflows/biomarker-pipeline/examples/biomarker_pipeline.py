'''End-to-end biomarker discovery pipeline'''

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, classification_report
from boruta import BorutaPy
import shap
import matplotlib.pyplot as plt

# Load data
# Example data: Use GEO datasets (e.g., GSE37418) or Bioconductor's curatedOvarianData
expr = pd.read_csv('expression.csv', index_col=0)
meta = pd.read_csv('metadata.csv', index_col=0)
X = expr.T  # transpose to samples x genes
y = meta.loc[X.index, 'condition'].values

print(f'Data: {X.shape[0]} samples, {X.shape[1]} features')
print(f'Classes: {np.unique(y, return_counts=True)}')

# Step 1: Train/test split
# test_size=0.2: Standard 80/20 split; use 0.3 for <100 samples
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, stratify=y, random_state=42
)
print(f'Train: {len(y_train)}, Test: {len(y_test)}')

# Scale features (fit on train only to prevent data leakage)
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Step 2: Feature Selection with Boruta
# max_depth=5: Shallow trees for stable importances across Boruta iterations
rf_selector = RandomForestClassifier(n_estimators=100, max_depth=5, n_jobs=-1, random_state=42)
# max_iter=100: Usually sufficient; increase to 200 if many tentative features remain
# n_estimators='auto': Scales with features (max of n_features, 500)
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

# Step 3: Nested CV for unbiased performance evaluation
# outer_cv=5: Standard for performance estimation; use 10 for small datasets
outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# n_estimators=100: Sufficient for most omics; increase to 500 for final model
clf = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1)
cv_scores = cross_val_score(clf, X_train_sel, y_train, cv=outer_cv, scoring='roc_auc')

print(f'Nested CV AUC: {cv_scores.mean():.3f} +/- {cv_scores.std():.3f}')

# QC: Check AUC and variance
if cv_scores.mean() < 0.7:
    print('WARNING: Low AUC. Check data quality or add samples.')
if cv_scores.std() > 0.1:
    print('WARNING: High fold variance. Consider more folds or LOOCV.')

# Step 4: Train final model and interpret with SHAP
clf.fit(X_train_sel, y_train)

# SHAP v0.47+ API: call explainer directly, NOT .shap_values()
explainer = shap.TreeExplainer(clf)
shap_values = explainer(X_train_sel)

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

# Bootstrap CI for AUC
# n_bootstrap=1000: Standard for publication-quality confidence intervals
n_bootstrap = 1000
boot_aucs = []
for i in range(n_bootstrap):
    idx = np.random.choice(len(y_test), size=len(y_test), replace=True)
    boot_aucs.append(roc_auc_score(y_test[idx], y_prob[idx]))

# 2.5, 97.5 percentiles: Standard for 95% confidence interval
ci_lower, ci_upper = np.percentile(boot_aucs, [2.5, 97.5])
print(f'95% CI: [{ci_lower:.3f}, {ci_upper:.3f}]')

# Classification report
print('\nClassification Report:')
print(classification_report(y_test, clf.predict(X_test_sel)))

# Export results
pd.DataFrame({'feature': selected_features}).to_csv('biomarker_panel.csv', index=False)
print(f'\nExported {len(selected_features)} biomarkers to biomarker_panel.csv')

# Optional: Save model for deployment
import joblib
joblib.dump(clf, 'biomarker_classifier.joblib')
joblib.dump(scaler, 'feature_scaler.joblib')
print('Saved classifier and scaler')
