import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import lombscargle, find_peaks
from statsmodels.stats.multitest import multipletests

np.random.seed(42)

# --- Simulate irregularly sampled periodic signal ---
# 30 observation times over 96h: irregular spacing mimics real experimental sampling
n_obs = 30
times = np.sort(np.random.uniform(0, 96, n_obs))

# True period = 18h: non-circadian oscillation (e.g., cell cycle in fast-dividing cells)
true_period = 18.0
amplitude = 2.0
mesor = 8.0
# SD = 0.5: moderate noise for normalized expression
signal = mesor + amplitude * np.sin(2 * np.pi * times / true_period) + np.random.normal(0, 0.5, n_obs)

signal_centered = signal - np.mean(signal)

# --- Lomb-Scargle periodogram ---
# Test periods from 6h to 72h
# min_period=6h: twice the median sampling interval (~3h) satisfies Nyquist
# max_period=72h: total duration / 1.3; need at least ~1.3 cycles for reliable detection
min_freq = 2 * np.pi / 72.0
max_freq = 2 * np.pi / 6.0
# 2000 frequency points: fine grid for precise period estimation
angular_freqs = np.linspace(min_freq, max_freq, 2000)
power = lombscargle(times, signal_centered, angular_freqs, normalize=True)

periods = 2 * np.pi / angular_freqs

# --- Find dominant peak ---
dominant_idx = np.argmax(power)
detected_period = periods[dominant_idx]
detected_power = power[dominant_idx]
print(f'True period: {true_period:.1f}h')
print(f'Detected period: {detected_period:.1f}h')
print(f'Peak power: {detected_power:.3f}')

# --- Significance via permutation ---
# n_perm=1000: standard for exploratory analysis; use 10000 for publication
n_perm = 1000
null_max_powers = np.zeros(n_perm)
for perm in range(n_perm):
    shuffled = np.random.permutation(signal_centered)
    null_power = lombscargle(times, shuffled, angular_freqs, normalize=True)
    null_max_powers[perm] = np.max(null_power)

# Empirical p-value: fraction of null peaks >= observed peak
empirical_p = np.mean(null_max_powers >= detected_power)
print(f'Empirical p-value (permutation): {empirical_p:.4f}')

# --- Significance levels from null distribution ---
# 99th percentile of null: FAP ~0.01 threshold
fap_01 = np.percentile(null_max_powers, 99)
# 95th percentile: FAP ~0.05 threshold
fap_05 = np.percentile(null_max_powers, 95)
print(f'FAP 0.01 power threshold: {fap_01:.3f}')
print(f'FAP 0.05 power threshold: {fap_05:.3f}')

# --- Multi-gene screening ---
n_genes = 200
n_periodic = 40

all_fap = []
all_periods = []
for g in range(n_genes):
    base = np.random.uniform(6, 12)
    if g < n_periodic:
        # Random period between 10-30h: simulates diverse oscillation periods
        gene_period = np.random.uniform(10, 30)
        gene_amp = np.random.uniform(1.0, 3.0)
        gene_signal = gene_amp * np.sin(2 * np.pi * times / gene_period) + np.random.normal(0, 0.5, n_obs)
    else:
        gene_signal = np.random.normal(0, 0.5, n_obs)

    gene_power = lombscargle(times, gene_signal, angular_freqs, normalize=True)
    max_power = np.max(gene_power)
    all_periods.append(periods[np.argmax(gene_power)])
    # Reuse null from single-gene permutation: valid when all genes share the same
    # observation times and frequency grid. For real data with per-gene missing values
    # or different sampling patterns, compute a separate null per gene.
    all_fap.append(np.mean(null_max_powers >= max_power))

# BH FDR correction for multiple testing
reject, qvals, _, _ = multipletests(all_fap, method='fdr_bh')
# q < 0.05: standard FDR threshold for periodic gene discovery
n_detected = np.sum(qvals < 0.05)
tp = np.sum(qvals[:n_periodic] < 0.05)
print(f'\nGenome-wide screening: {n_detected}/{n_genes} periodic genes detected')
print(f'True positives: {tp}/{n_periodic}')

# --- Visualization ---
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

axes[0, 0].scatter(times, signal, c='steelblue', s=40, zorder=3, label='Observed')
t_fine = np.linspace(0, 96, 500)
axes[0, 0].plot(t_fine, mesor + amplitude * np.sin(2 * np.pi * t_fine / true_period),
                'r--', alpha=0.5, label=f'True ({true_period}h)')
axes[0, 0].set_xlabel('Time (hours)')
axes[0, 0].set_ylabel('Expression')
axes[0, 0].set_title('Irregularly sampled periodic signal')
axes[0, 0].legend()

axes[0, 1].plot(periods, power, 'steelblue', linewidth=1.5)
axes[0, 1].axvline(true_period, color='red', linestyle='--', alpha=0.7, label=f'True period ({true_period}h)')
axes[0, 1].axhline(fap_01, color='gray', linestyle=':', alpha=0.7, label='FAP 0.01')
axes[0, 1].axhline(fap_05, color='gray', linestyle='--', alpha=0.5, label='FAP 0.05')
axes[0, 1].set_xlabel('Period (hours)')
axes[0, 1].set_ylabel('Normalized power')
axes[0, 1].set_title(f'Lomb-Scargle periodogram (peak={detected_period:.1f}h)')
axes[0, 1].legend(fontsize=8)
axes[0, 1].invert_xaxis()

axes[1, 0].hist(null_max_powers, bins=30, color='gray', edgecolor='black', alpha=0.7, density=True, label='Null')
axes[1, 0].axvline(detected_power, color='red', linewidth=2, label=f'Observed (p={empirical_p:.3f})')
axes[1, 0].set_xlabel('Maximum power')
axes[1, 0].set_ylabel('Density')
axes[1, 0].set_title('Permutation null distribution')
axes[1, 0].legend()

periodic_mask = qvals < 0.05
axes[1, 1].scatter(np.array(all_periods)[~periodic_mask], -np.log10(np.array(qvals)[~periodic_mask]),
                   c='gray', s=20, alpha=0.5, label='Non-significant')
axes[1, 1].scatter(np.array(all_periods)[periodic_mask], -np.log10(np.array(qvals)[periodic_mask]),
                   c='coral', s=30, alpha=0.8, label='Significant (q<0.05)')
axes[1, 1].axhline(-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
axes[1, 1].set_xlabel('Dominant period (hours)')
axes[1, 1].set_ylabel('-log10(q-value)')
axes[1, 1].set_title(f'Genome-wide periodicity ({n_detected} detected)')
axes[1, 1].legend(fontsize=8)

plt.tight_layout()
plt.savefig('lomb_scargle_results.png', dpi=150, bbox_inches='tight')
plt.show()
print('\nPlot saved to lomb_scargle_results.png')
