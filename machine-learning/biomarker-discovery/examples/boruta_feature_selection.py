'''Boruta all-relevant selection keeps the whole correlated module (unlike LASSO).

Runs end-to-end on synthetic data. Five co-expressed genes (g0-g4) all carry the
signal. An all-relevant selector (Boruta) confirms the redundant module members,
whereas a minimal-optimal selector would keep only one -- the distinction that
decides whether absence from a list means "irrelevant" or merely "redundant".
'''
# Reference: numpy 1.26+, pandas 2.2+, scikit-learn 1.4+, boruta 0.4+ | Verify API if version differs

import numpy as np
import pandas as pd
from boruta import BorutaPy
from sklearn.ensemble import RandomForestClassifier

rng = np.random.default_rng(0)
n, p = 200, 60
latent = rng.normal(size=n)                            # shared biological signal
X = rng.normal(size=(n, p))
X[:, :5] = latent[:, None] + rng.normal(scale=0.3, size=(n, 5))   # g0-g4: one co-expressed module
y = (latent + rng.normal(scale=0.3, size=n) > 0).astype(int)
X = pd.DataFrame(X, columns=[f'g{i}' for i in range(p)])

# Boruta needs a tree estimator and numpy arrays. perc=100 uses the max shadow importance.
rf = RandomForestClassifier(n_estimators=200, n_jobs=-1, max_depth=5, random_state=42)
boruta = BorutaPy(rf, n_estimators='auto', perc=100, two_step=True, max_iter=100, random_state=42)
boruta.fit(X.values, y)                                # numpy arrays, not pandas

confirmed = list(X.columns[boruta.support_])
tentative = list(X.columns[boruta.support_weak_])
print(f'Confirmed all-relevant ({len(confirmed)}): {confirmed}')
print(f'Tentative ({len(tentative)}): {tentative}')
module_recovered = sum(g in confirmed for g in ['g0', 'g1', 'g2', 'g3', 'g4'])
print(f'\nModule members g0-g4 confirmed: {module_recovered}/5 '
      f'(all-relevant keeps the redundant module; a minimal-optimal selector would keep ~1)')
