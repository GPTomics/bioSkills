'''SHAP interpretation for omics classifiers (v0.47+ API)'''

import pandas as pd
import numpy as np
import shap
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

expr = pd.read_csv('expression.csv', index_col=0)
meta = pd.read_csv('metadata.csv', index_col=0)

X = expr.T
y = meta.loc[X.index, 'condition'].values

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, stratify=y, random_state=42)

model = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1)
model.fit(X_train, y_train)

# SHAP v0.47+ API: call explainer directly, NOT .shap_values()
explainer = shap.TreeExplainer(model)
shap_values = explainer(X_test)
print(f'SHAP values shape: {shap_values.values.shape}')

# Beeswarm: shows importance AND direction of effect
# max_display=20: Top 20 features; increase for more
shap.plots.beeswarm(shap_values, max_display=20, show=False)
plt.tight_layout()
plt.savefig('shap_beeswarm.png', dpi=150, bbox_inches='tight')
plt.close()

# Bar plot: mean |SHAP|
shap.plots.bar(shap_values, max_display=20, show=False)
plt.tight_layout()
plt.savefig('shap_bar.png', dpi=150, bbox_inches='tight')
plt.close()

# Waterfall for single prediction
sample_idx = 0
shap.plots.waterfall(shap_values[sample_idx], max_display=15, show=False)
plt.tight_layout()
plt.savefig('shap_waterfall_sample0.png', dpi=150, bbox_inches='tight')
plt.close()

# Extract feature importances
mean_shap = np.abs(shap_values.values).mean(axis=0)
feature_importance = pd.DataFrame({
    'feature': X_test.columns,
    'mean_abs_shap': mean_shap
}).sort_values('mean_abs_shap', ascending=False)

feature_importance.to_csv('shap_feature_importance.csv', index=False)
print('\nTop 20 features by SHAP importance:')
feature_importance.head(20)
