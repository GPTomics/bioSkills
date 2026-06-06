'''LIME is non-reproducible: the same prediction yields different top features per seed.

Runs end-to-end on synthetic data. Explaining ONE fixed prediction with LIME under
several random seeds returns different top-feature sets, so a LIME ranking is not a
stable finding. Use LIME only to eyeball one local prediction, never for global ranking.
'''
# Reference: numpy 1.26+, scikit-learn 1.4+, lime 0.2+ | Verify API if version differs

import numpy as np
from lime.lime_tabular import LimeTabularExplainer
from sklearn.datasets import make_classification
from sklearn.ensemble import RandomForestClassifier

X, y = make_classification(n_samples=400, n_features=30, n_informative=10, random_state=0)
model = RandomForestClassifier(n_estimators=200, random_state=0).fit(X, y)
instance = X[0]                                          # one fixed prediction to explain

top_sets = []
for seed in range(5):
    # Same instance, same model -- only the LIME random seed changes.
    explainer = LimeTabularExplainer(X, mode='classification',
                                     discretize_continuous=True, random_state=seed)
    exp = explainer.explain_instance(instance, model.predict_proba, num_features=5, num_samples=5000)
    top_feats = tuple(sorted(idx for idx, _ in exp.local_exp[exp.available_labels()[0]][:5]))
    top_sets.append(top_feats)
    print(f'seed {seed}: top-5 features = {top_feats}')

distinct = len(set(top_sets))
print(f'\n{distinct} distinct top-5 sets across 5 seeds for the SAME prediction.')
print('LIME rankings are seed-dependent -- never use them as a stable importance measure.')
