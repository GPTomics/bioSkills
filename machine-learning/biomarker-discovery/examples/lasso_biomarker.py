'''LASSO feature selection with stability analysis'''

import pandas as pd
import numpy as np
from sklearn.linear_model import LassoCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold

expr = pd.read_csv('expression.csv', index_col=0)
meta = pd.read_csv('metadata.csv', index_col=0)

X = expr.T
y = meta.loc[X.index, 'condition'].values

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# cv=5: Standard 5-fold CV balances bias-variance; use 10 for smaller datasets
lasso = LassoCV(cv=5, random_state=42, n_jobs=-1)
lasso.fit(X_scaled, y)

selected_mask = lasso.coef_ != 0
print(f'LASSO selected {selected_mask.sum()} features at alpha={lasso.alpha_:.4f}')

coefs = pd.DataFrame({
    'feature': X.columns,
    'coefficient': lasso.coef_,
    'selected': selected_mask
})
coefs = coefs.sort_values('coefficient', key=abs, ascending=False)
coefs.to_csv('lasso_features.csv', index=False)

# Stability selection via bootstrap
# n_bootstrap=100: Sufficient for stable estimates; 500-1000 for publication-quality
n_bootstrap = 100
selection_counts = np.zeros(X.shape[1])

print(f'\nRunning {n_bootstrap} bootstrap iterations...')
for i in range(n_bootstrap):
    idx = np.random.choice(len(X), size=len(X), replace=True)
    X_boot = X_scaled[idx]
    y_boot = y[idx]

    # cv=3: Fewer folds for speed in bootstrap; repeated sampling provides stability
    lasso_boot = LassoCV(cv=3, random_state=i, n_jobs=-1)
    lasso_boot.fit(X_boot, y_boot)
    selection_counts += (lasso_boot.coef_ != 0)

# threshold=0.6: Feature selected in >60% of bootstraps is stable
# Lower (0.5) for more features, higher (0.8) for stricter selection
stability = selection_counts / n_bootstrap
stable_mask = stability > 0.6

stability_df = pd.DataFrame({
    'feature': X.columns,
    'selection_frequency': stability,
    'stable': stable_mask
}).sort_values('selection_frequency', ascending=False)

print(f'Stable features (>60% selection): {stable_mask.sum()}')
stability_df.to_csv('lasso_stability.csv', index=False)
stability_df[stability_df['stable']].head(20)
