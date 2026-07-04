# Reference: numpy 1.26+, ruptures 1.1+, matplotlib 3.8+ | Verify API if version differs
import os
import tempfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import ruptures as rpt

np.random.seed(42)

# --- Simulate a piecewise expression signal ---
n_timepoints = 20
timepoints = np.linspace(0, 48, n_timepoints)
# 3 regimes: basal (0-16h), induced (16-32h, step up), recovery (32-48h, partial return).
signal_clean = np.piecewise(
    timepoints,
    [timepoints < 16, (timepoints >= 16) & (timepoints < 32), timepoints >= 32],
    [lambda t: 8.0 + 0.05 * t,
     lambda t: 8.0 + 0.05 * 16 + 2.0,
     lambda t: 8.0 + 0.05 * 16 + 0.5])
# SD=0.3: moderate noise on a log-expression scale.
signal = signal_clean + np.random.normal(0, 0.3, n_timepoints)

# --- Pelt changepoint detection ---
# model='l2' detects changes in MEAN (piecewise-constant level), the natural cost for step-like
# regimes and the cost for which the BIC penalty log(n)*sigma^2 is derived (rbf would mismatch it).
# min_size=2: minimum segment length. Pelt is exact; the penalty IS the number-of-breaks knob.
n = len(signal)
# Noise variance from lag-1 differences, NOT total np.var(signal): the total includes
# between-regime variation, which inflates the penalty and UNDER-detects real changepoints.
sigma2 = np.var(np.diff(signal)) / 2.0
penalty_bic = np.log(n) * sigma2
bkps_pelt = rpt.Pelt(model='l2', min_size=2).fit(signal).predict(pen=penalty_bic)
# predict returns break indices INCLUDING the terminal index n.
print(f'Pelt changepoints (BIC penalty={penalty_bic:.3f}): {bkps_pelt}')
print(f'Number of changepoints: {len(bkps_pelt) - 1}')

# --- Binary segmentation (approximate, needs a known number of breaks) ---
bkps_binseg = rpt.Binseg(model='l2', min_size=2).fit(signal).predict(n_bkps=2)
print(f'BinSeg changepoints (n_bkps=2): {bkps_binseg}')

# --- Sensitivity to the penalty: the penalty choice IS the number-of-changepoints choice ---
algo_pelt = rpt.Pelt(model='l2', min_size=2).fit(signal)
for pen in [0.1, 0.3, 0.5, 1.0, 2.0]:
    bkps = algo_pelt.predict(pen=pen)
    print(f'  Penalty {pen:4.1f}: {len(bkps) - 1} changepoints at indices {bkps[:-1]}')

# --- Multi-gene changepoint detection ---
n_genes, n_with_changepoints = 100, 30
all_signals = np.zeros((n_genes, n_timepoints))
for i in range(n_genes):
    base = np.random.uniform(6, 12)
    if i < n_with_changepoints:
        cp_time = np.random.choice(range(5, 15))
        all_signals[i, :cp_time] = base
        all_signals[i, cp_time:] = base + np.random.uniform(1.0, 3.0)
    else:
        all_signals[i, :] = base
    all_signals[i, :] += np.random.normal(0, 0.3, n_timepoints)

gene_results = []
for i in range(n_genes):
    sig = all_signals[i]
    pen = np.log(n_timepoints) * (np.var(np.diff(sig)) / 2.0)
    bkps = rpt.Pelt(model='l2', min_size=2).fit(sig).predict(pen=pen)
    gene_results.append({'gene': f'gene_{i}', 'n_changepoints': len(bkps) - 1,
                         'changepoints': bkps[:-1]})

genes_with_changes = [r for r in gene_results if r['n_changepoints'] > 0]
print(f'\nGenome-wide: {len(genes_with_changes)}/{n_genes} genes with detected changepoints')
print(f'True positives in first {n_with_changepoints}: '
      f'{sum(1 for r in gene_results[:n_with_changepoints] if r["n_changepoints"] > 0)}')

# --- Visualization (written to a temp file, then removed: no stray outputs) ---
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
axes[0, 0].scatter(timepoints, signal, c='steelblue', s=40, zorder=3)
axes[0, 0].plot(timepoints, signal_clean, 'k--', alpha=0.5, label='True signal')
for bkp in bkps_pelt[:-1]:
    axes[0, 0].axvline(timepoints[bkp - 1], color='red', linestyle='--', linewidth=1.5, alpha=0.8)
axes[0, 0].set(xlabel='Time (hours)', ylabel='Expression',
               title=f'Pelt: {len(bkps_pelt) - 1} changepoints detected')
axes[0, 0].legend()

axes[0, 1].scatter(timepoints, signal, c='steelblue', s=40, zorder=3)
for bkp in bkps_binseg[:-1]:
    axes[0, 1].axvline(timepoints[bkp - 1], color='orange', linestyle='--', linewidth=1.5, alpha=0.8)
axes[0, 1].set(xlabel='Time (hours)', ylabel='Expression',
               title=f'BinSeg: {len(bkps_binseg) - 1} changepoints detected')

n_changes_dist = [r['n_changepoints'] for r in gene_results]
axes[1, 0].hist(n_changes_dist, bins=range(max(n_changes_dist) + 2), color='coral',
                edgecolor='black', align='left')
axes[1, 0].set(xlabel='Number of changepoints', ylabel='Number of genes',
               title='Changepoint count distribution (genome-wide)')

axes[1, 1].scatter(timepoints, all_signals[0], c='steelblue', s=40, zorder=3)
for bkp in gene_results[0]['changepoints']:
    axes[1, 1].axvline(timepoints[bkp - 1], color='red', linestyle='--', linewidth=1.5)
axes[1, 1].set(xlabel='Time (hours)', ylabel='Expression',
               title=f"gene_0: {len(gene_results[0]['changepoints'])} changepoint(s)")

plt.tight_layout()
out_png = tempfile.NamedTemporaryFile(suffix='.png', delete=False).name
plt.savefig(out_png, dpi=150, bbox_inches='tight')
plt.close(fig)
os.unlink(out_png)
print('\nDone (figure rendered to a temp file and removed).')
