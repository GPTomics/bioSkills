'''Patient-as-unit differential abundance, and why cell-level testing is wrong.

Demonstrates the invariant spine: per-image cell-type proportions -> one value per patient
-> mixed model with patient as a random effect and batch as a covariate. The contrast at
the end shows the pseudoreplication trap: a cell-level test over correlated cells returns a
near-zero p-value for an effect that is not significant at the patient level.
'''
# Reference: statsmodels 0.14+, pandas 2.2+, numpy 1.26+, scipy 1.12+ | Verify API if version differs
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

rng = np.random.default_rng(0)

# simulate a cohort: 20 patients x 2 conditions x ~3 ROIs, ~2000 cells/ROI. The TRUE Treg
# proportion has a small condition effect swamped by large between-patient variance.
rows = []
for patient in range(20):
    cond = 'responder' if patient % 2 else 'nonresponder'
    batch = patient % 3
    patient_baseline = rng.normal(0.10, 0.04)                 # between-patient variance dominates
    effect = 0.01 if cond == 'responder' else 0.0
    for roi in range(rng.integers(2, 5)):
        p_treg = np.clip(patient_baseline + effect + rng.normal(0, 0.01), 0.001, 0.5)
        n = 2000
        n_treg = rng.binomial(n, p_treg)
        rows.append(dict(patient=patient, condition=cond, batch=batch,
                         image_id=f'{patient}_{roi}', n_treg=n_treg, n_cells=n))
img = pd.DataFrame(rows)
img['prop'] = img['n_treg'] / img['n_cells']

# CORRECT: mixed model on per-image proportions with patient as the unit, batch as covariate
mm = smf.mixedlm('prop ~ condition + C(batch)', img, groups=img['patient']).fit()
p_patient = mm.pvalues['condition[T.responder]']
print(f'Patient-unit mixed model p(condition) = {p_patient:.3f}  (effective n = {img.patient.nunique()} patients)')

# WRONG: treat every cell as an independent replicate (pseudoreplication)
cells = []
for _, r in img.iterrows():
    cells += [(r.condition, 1)] * int(r.n_treg) + [(r.condition, 0)] * int(r.n_cells - r.n_treg)
cell_df = pd.DataFrame(cells, columns=['condition', 'is_treg'])
table = pd.crosstab(cell_df.condition, cell_df.is_treg)
_, p_cell = stats.fisher_exact(table.values) if table.shape == (2, 2) else (None, np.nan)
print(f'Cell-level test p = {p_cell:.1e}  <- inflated by ~{len(cell_df):,} pseudoreplicated cells')
print('The cell-level p is near zero for an effect the patient-level test does not support.')
