import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cosinorpy import cosinor, cosinor1, file_parser
from statsmodels.stats.multitest import multipletests

np.random.seed(42)

# --- Simulate circadian expression data ---
# 48h sampled every 4h = 13 timepoints covering 2 full 24h cycles
timepoints = np.arange(0, 48 + 1, 4)
n_genes = 200
n_rhythmic = 50

expression_records = []
for i in range(n_genes):
    mesor = np.random.uniform(5, 12)
    if i < n_rhythmic:
        # 24h period: standard circadian oscillation
        amplitude = np.random.uniform(1.0, 3.0)
        # Phase drawn uniformly across 24h cycle
        phase = np.random.uniform(0, 2 * np.pi)
        values = mesor + amplitude * np.cos(2 * np.pi * timepoints / 24 - phase)
    else:
        amplitude = 0
        values = np.full(len(timepoints), mesor)
    # Gaussian noise with SD ~0.5: typical for normalized RNA-seq log-expression
    noise = np.random.normal(0, 0.5, len(timepoints))
    for t_idx, t in enumerate(timepoints):
        expression_records.append({
            'test': f'gene_{i}',
            'x': t,
            'y': values[t_idx] + noise[t_idx]
        })

df = pd.DataFrame(expression_records)

# --- Single-gene cosinor fit ---
# period=24: testing for standard 24h circadian rhythm
# fit_group works on a DataFrame and returns a DataFrame (simpler than fit_me for batch use)
single_result = cosinor.fit_group(df[df['test'] == 'gene_0'], period=24, n_components=1)
print('Single gene cosinor result:')
print(single_result[['test', 'amplitude', 'acrophase', 'p']])

# --- Genome-wide cosinor analysis ---
# fit_group handles all genes at once when given multi-gene long-format DataFrame
results_df = cosinor.fit_group(df, period=24, n_components=1)

# --- Multiple testing correction ---
# BH FDR: standard correction for genome-wide rhythmicity testing
valid_mask = results_df['p'].notna()
reject, qvals, _, _ = multipletests(results_df.loc[valid_mask, 'p'], method='fdr_bh')
results_df.loc[valid_mask, 'q_bh'] = qvals

# q < 0.05: standard FDR threshold for circadian gene discovery
rhythmic_genes = results_df[results_df['q_bh'] < 0.05].copy()
print(f'\nRhythmic genes detected: {len(rhythmic_genes)} / {n_genes}')
print(f'True positives in first {n_rhythmic} genes: {sum(rhythmic_genes["test"].str.extract(r"(\d+)")[0].astype(int) < n_rhythmic)}')

# --- Population-mean cosinor (multi-subject) ---
# Simulate 3 subjects with shared rhythm parameters
subject_records = []
for subj in range(3):
    mesor = 8.0 + np.random.normal(0, 0.3)
    # Amplitude ~2.0: moderate circadian amplitude for demonstration
    amplitude = 2.0 + np.random.normal(0, 0.2)
    phase = 1.0 + np.random.normal(0, 0.1)
    values = mesor + amplitude * np.cos(2 * np.pi * timepoints / 24 - phase)
    noise = np.random.normal(0, 0.4, len(timepoints))
    for t_idx, t in enumerate(timepoints):
        subject_records.append({
            'test': f'subject_{subj}',
            'x': t,
            'y': values[t_idx] + noise[t_idx]
        })

pop_df = pd.DataFrame(subject_records)
pop_result = cosinor1.population_fit_cosinor(pop_df, period=24)
print('\nPopulation-mean cosinor:')
print(pop_result)

# --- Visualization ---
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

gene_data = df[df['test'] == 'gene_0']
axes[0].scatter(gene_data['x'], gene_data['y'], c='steelblue', s=40, zorder=3)
t_fine = np.linspace(0, 48, 200)
row = results_df[results_df['test'] == 'gene_0'].iloc[0]
fitted = row['mesor'] + row['amplitude'] * np.cos(2 * np.pi * t_fine / 24 - row['acrophase'])
axes[0].plot(t_fine, fitted, 'r-', linewidth=2)
axes[0].set_xlabel('Time (hours)')
axes[0].set_ylabel('Expression')
axes[0].set_title(f'gene_0 cosinor fit (p={row["p"]:.2e})')

if len(rhythmic_genes) > 0:
    acrophases_hours = (rhythmic_genes['acrophase'] * 24 / (2 * np.pi)) % 24
    axes[1].hist(acrophases_hours, bins=24, range=(0, 24), color='coral', edgecolor='black')
    axes[1].set_xlabel('Acrophase (hours)')
    axes[1].set_ylabel('Number of genes')
    axes[1].set_title('Phase distribution of rhythmic genes')

plt.tight_layout()
plt.savefig('cosinor_results.png', dpi=150, bbox_inches='tight')
plt.show()
print('\nPlot saved to cosinor_results.png')
