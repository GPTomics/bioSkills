# Reference: CosinorPy 3.1+ (requires numpy<2.0), pandas 2.2+, statsmodels 0.14+ | Verify API if version differs
# CosinorPy 3.1 imports as `from CosinorPy import ...` (capitalized package); 0.x/1.x used lowercase `cosinorpy`.
import os
import tempfile
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from CosinorPy import cosinor, cosinor1
from statsmodels.stats.multitest import multipletests

np.random.seed(42)

# --- Simulate circadian expression data ---
# 48h sampled every 4h = 13 timepoints covering 2 full 24h cycles (>=2 cycles, >=6/cycle design minima).
timepoints = np.arange(0, 48 + 1, 4)
n_genes = 200
n_rhythmic = 50

records = []
for i in range(n_genes):
    mesor = np.random.uniform(5, 12)
    if i < n_rhythmic:
        amplitude = np.random.uniform(1.0, 3.0)          # moderate circadian amplitude
        phase = np.random.uniform(0, 2 * np.pi)          # phases drawn uniformly across the cycle
        values = mesor + amplitude * np.cos(2 * np.pi * timepoints / 24 - phase)
    else:
        values = np.full(len(timepoints), mesor)
    noise = np.random.normal(0, 0.5, len(timepoints))    # SD 0.5: typical normalized log-expression noise
    for t, y in zip(timepoints, values + noise):
        records.append({'test': f'gene_{i}', 'x': t, 'y': y})

df = pd.DataFrame(records)

# --- Genome-wide cosinor with built-in q-values ---
# fit_group returns a DataFrame whose 'q' column is ALREADY BH-adjusted across the fitted group.
results = cosinor.fit_group(df, period=24, n_components=1, plot=False)

# Recompute BH only to control the correction SET explicitly; default multipletests method is Holm-Sidak,
# so fdr_bh must be passed for Benjamini-Hochberg.
valid = results['p'].notna()
results.loc[valid, 'q_bh'] = multipletests(results.loc[valid, 'p'], method='fdr_bh')[1]

# rAMP = amplitude / MESOR: comparable across genes because it normalizes out expression level.
# rAMP > 0.1 (>=10% of baseline): conventional biological-relevance floor - sweep it, do not treat as law.
results['rAMP'] = results['amplitude'] / results['mesor']
rhythmic = results[(results['q'] < 0.05) & (results['rAMP'] > 0.1)].copy()

print('fit_group columns:', list(results.columns))
true_pos = (rhythmic['test'].str.extract(r'(\d+)')[0].astype(int) < n_rhythmic).sum()
print(f'Rhythmic (q<0.05 & rAMP>0.1): {len(rhythmic)} / {n_genes}; true positives: {true_pos} / {n_rhythmic}')

# --- Population-mean cosinor (multi-subject) ---
# Three subjects sharing rhythm parameters; population_fit_cosinor returns a DICT with group CIs and p-values,
# correctly propagating BETWEEN-subject variance rather than pooling all points into one fit.
subject_records = []
for subj in range(3):
    mesor = 8.0 + np.random.normal(0, 0.3)
    amplitude = 2.0 + np.random.normal(0, 0.2)
    phase = 1.0 + np.random.normal(0, 0.1)
    values = mesor + amplitude * np.cos(2 * np.pi * timepoints / 24 - phase)
    noise = np.random.normal(0, 0.4, len(timepoints))
    for t, y in zip(timepoints, values + noise):
        subject_records.append({'test': f'subject_{subj}', 'x': t, 'y': y})

pop = cosinor1.population_fit_cosinor(pd.DataFrame(subject_records), period=24, plot_on=False)
print('\nPopulation-mean cosinor: group amp CI', pop['confint']['amp'], 'p_amp', round(pop['p_amp'], 4))

# --- Visualization (written to a temp dir and cleaned up so no stray files remain) ---
# CosinorPy stores acrophase as phi in y = M + A*cos(2*pi*t/T + phi) = atan2(-gamma, beta) (usually negative),
# so the fitted curve uses +acrophase and peak-hour = (-acrophase)*T/(2*pi) mod T.
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
row = results[results['test'] == 'gene_0'].iloc[0]
g0 = df[df['test'] == 'gene_0']
t_fine = np.linspace(0, 48, 200)
fitted = row['mesor'] + row['amplitude'] * np.cos(2 * np.pi * t_fine / 24 + row['acrophase'])
axes[0].scatter(g0['x'], g0['y'], c='steelblue', s=40, zorder=3)
axes[0].plot(t_fine, fitted, 'r-', linewidth=2)
axes[0].set(xlabel='Time (hours)', ylabel='Expression', title=f'gene_0 cosinor fit (q={row["q"]:.1e})')

peak_hours = ((-rhythmic['acrophase']) * 24 / (2 * np.pi)) % 24
axes[1].hist(peak_hours, bins=24, range=(0, 24), color='coral', edgecolor='black')
axes[1].set(xlabel='Peak time (hours)', ylabel='Number of genes', title='Phase distribution of rhythmic genes')

plt.tight_layout()
out_path = os.path.join(tempfile.mkdtemp(), 'cosinor_results.png')
plt.savefig(out_path, dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'\nPlot written to {out_path}')
os.remove(out_path)
os.rmdir(os.path.dirname(out_path))
