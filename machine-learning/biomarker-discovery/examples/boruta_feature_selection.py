'''Boruta all-relevant feature selection for biomarker discovery'''

import pandas as pd
import numpy as np
from boruta import BorutaPy
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, f_classif

expr = pd.read_csv('expression.csv', index_col=0)
meta = pd.read_csv('metadata.csv', index_col=0)

X = expr.T
y = meta.loc[X.index, 'condition'].values
print(f'Data: {X.shape[0]} samples, {X.shape[1]} features')

# Pre-filter to top 5000 by ANOVA F-statistic
# k=5000: Reduces computation; adjust based on total features
if X.shape[1] > 5000:
    selector = SelectKBest(f_classif, k=5000)
    selector.fit(X, y)
    X_filtered = X.iloc[:, selector.get_support()]
    print(f'Pre-filtered to {X_filtered.shape[1]} features')
else:
    X_filtered = X

# max_depth=5: Shallow trees for speed and stable feature importances across Boruta iterations
rf = RandomForestClassifier(n_estimators=100, n_jobs=-1, max_depth=5, random_state=42)

# max_iter=100: Usually sufficient; set to 200 if many tentative features remain
# n_estimators='auto': Scales with features (max of n_features, 500)
boruta = BorutaPy(rf, n_estimators='auto', max_iter=100, random_state=42, verbose=2)
boruta.fit(X_filtered.values, y)

results = pd.DataFrame({
    'feature': X_filtered.columns,
    'rank': boruta.ranking_,
    'selected': boruta.support_,
    'tentative': boruta.support_weak_
}).sort_values('rank')

selected = results[results['selected']]['feature'].tolist()
tentative = results[results['tentative']]['feature'].tolist()

print(f'\nSelected: {len(selected)} features')
print(f'Tentative: {len(tentative)} features')
print(f'Rejected: {len(results) - len(selected) - len(tentative)} features')

results.to_csv('boruta_results.csv', index=False)
print('\nTop 20 selected features:')
results[results['selected']].head(20)
