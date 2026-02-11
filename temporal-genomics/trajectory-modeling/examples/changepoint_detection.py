import numpy as np
import matplotlib.pyplot as plt
import ruptures as rpt

np.random.seed(42)

# --- Simulate piecewise expression signal ---
# 20 timepoints: dense enough for changepoint detection
n_timepoints = 20
timepoints = np.linspace(0, 48, n_timepoints)

# 3 regimes: basal (0-16h), induced (16-32h), recovery (32-48h)
signal_clean = np.piecewise(
    timepoints,
    [timepoints < 16, (timepoints >= 16) & (timepoints < 32), timepoints >= 32],
    [lambda t: 8.0 + 0.05 * t,          # basal: slow rise
     lambda t: 8.0 + 0.05 * 16 + 2.0,   # induced: step up by 2
     lambda t: 8.0 + 0.05 * 16 + 0.5]   # recovery: partial return
)
# SD = 0.3: moderate noise for log-expression
signal = signal_clean + np.random.normal(0, 0.3, n_timepoints)

# --- Pelt changepoint detection ---
# model='rbf': radial basis function kernel; detects changes in mean and variance
# min_size=2: minimum 2 timepoints per segment; prevents trivially small segments
algo_pelt = rpt.Pelt(model='rbf', min_size=2).fit(signal)

# BIC penalty: log(n) * variance; standard model selection criterion
# Balances goodness-of-fit against model complexity (number of changepoints)
n = len(signal)
penalty_bic = np.log(n) * np.var(signal)
bkps_pelt = algo_pelt.predict(pen=penalty_bic)
print(f'Pelt changepoints (BIC penalty={penalty_bic:.3f}): {bkps_pelt}')
print(f'Number of changepoints: {len(bkps_pelt) - 1}')

# --- Binary segmentation (faster alternative) ---
# n_bkps=2: expected number of changepoints based on experimental design
# BinSeg is O(n log n) vs Pelt O(n); use for very long series
algo_binseg = rpt.Binseg(model='rbf', min_size=2).fit(signal)
bkps_binseg = algo_binseg.predict(n_bkps=2)
print(f'BinSeg changepoints (n_bkps=2): {bkps_binseg}')

# --- Sensitivity to penalty ---
# Lower penalty = more changepoints (liberal); higher = fewer (conservative)
penalties = [0.5, 1.0, 2.0, 5.0, 10.0]
for pen in penalties:
    bkps = algo_pelt.predict(pen=pen)
    print(f'  Penalty {pen:5.1f}: {len(bkps) - 1} changepoints at indices {bkps[:-1]}')

# --- Multi-gene changepoint detection ---
n_genes = 100
n_with_changepoints = 30

all_signals = np.zeros((n_genes, n_timepoints))
for i in range(n_genes):
    base = np.random.uniform(6, 12)
    if i < n_with_changepoints:
        cp_time = np.random.choice(range(5, 15))
        shift = np.random.uniform(1.0, 3.0)
        all_signals[i, :cp_time] = base
        all_signals[i, cp_time:] = base + shift
    else:
        all_signals[i, :] = base
    all_signals[i, :] += np.random.normal(0, 0.3, n_timepoints)

gene_results = []
for i in range(n_genes):
    algo = rpt.Pelt(model='rbf', min_size=2).fit(all_signals[i])
    pen = np.log(n_timepoints) * np.var(all_signals[i])
    bkps = algo.predict(pen=pen)
    n_bkps = len(bkps) - 1
    gene_results.append({'gene': f'gene_{i}', 'n_changepoints': n_bkps,
                          'changepoints': bkps[:-1]})

genes_with_changes = [r for r in gene_results if r['n_changepoints'] > 0]
print(f'\nGenome-wide: {len(genes_with_changes)}/{n_genes} genes with detected changepoints')
print(f'True positives in first {n_with_changepoints}: '
      f'{sum(1 for r in gene_results[:n_with_changepoints] if r["n_changepoints"] > 0)}')

# --- Visualization ---
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

axes[0, 0].scatter(timepoints, signal, c='steelblue', s=40, zorder=3)
axes[0, 0].plot(timepoints, signal_clean, 'k--', alpha=0.5, label='True signal')
for bkp in bkps_pelt[:-1]:
    axes[0, 0].axvline(timepoints[bkp - 1], color='red', linestyle='--', linewidth=1.5, alpha=0.8)
axes[0, 0].set_xlabel('Time (hours)')
axes[0, 0].set_ylabel('Expression')
axes[0, 0].set_title(f'Pelt: {len(bkps_pelt) - 1} changepoints detected')
axes[0, 0].legend()

axes[0, 1].scatter(timepoints, signal, c='steelblue', s=40, zorder=3)
for bkp in bkps_binseg[:-1]:
    axes[0, 1].axvline(timepoints[bkp - 1], color='orange', linestyle='--', linewidth=1.5, alpha=0.8)
axes[0, 1].set_xlabel('Time (hours)')
axes[0, 1].set_ylabel('Expression')
axes[0, 1].set_title(f'BinSeg: {len(bkps_binseg) - 1} changepoints detected')

n_changes_dist = [r['n_changepoints'] for r in gene_results]
axes[1, 0].hist(n_changes_dist, bins=range(max(n_changes_dist) + 2),
                color='coral', edgecolor='black', align='left')
axes[1, 0].set_xlabel('Number of changepoints')
axes[1, 0].set_ylabel('Number of genes')
axes[1, 0].set_title('Changepoint count distribution (genome-wide)')

example_idx = 0
axes[1, 1].scatter(timepoints, all_signals[example_idx], c='steelblue', s=40, zorder=3)
example_bkps = gene_results[example_idx]['changepoints']
for bkp in example_bkps:
    axes[1, 1].axvline(timepoints[bkp - 1], color='red', linestyle='--', linewidth=1.5)
axes[1, 1].set_xlabel('Time (hours)')
axes[1, 1].set_ylabel('Expression')
axes[1, 1].set_title(f'gene_0: {len(example_bkps)} changepoint(s)')

plt.tight_layout()
plt.savefig('changepoint_results.png', dpi=150, bbox_inches='tight')
plt.show()
print('\nPlot saved to changepoint_results.png')
