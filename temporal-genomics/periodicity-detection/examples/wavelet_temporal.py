import numpy as np
import pywt
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

np.random.seed(42)

# --- Simulate signal with time-varying periodicity ---
# 200 evenly spaced timepoints over 96h (sampling interval ~0.48h)
n_timepoints = 200
times = np.linspace(0, 96, n_timepoints)
sampling_interval = times[1] - times[0]

# Phase 1 (0-48h): 12h oscillation (ultradian rhythm)
# Phase 2 (48-96h): 24h oscillation (circadian rhythm)
# Transition simulates entrainment shift or developmental change
signal = np.zeros(n_timepoints)
for i, t in enumerate(times):
    if t < 48:
        # 12h period with amplitude 2.0
        signal[i] = 2.0 * np.sin(2 * np.pi * t / 12.0)
    else:
        # 24h period with amplitude 1.5 (slightly lower; common after entrainment shift)
        signal[i] = 1.5 * np.sin(2 * np.pi * t / 24.0)

# SD = 0.4: moderate noise for time-frequency analysis
noise = np.random.normal(0, 0.4, n_timepoints)
signal_noisy = signal + noise

# --- Continuous Wavelet Transform ---
# cmor1.5-1.0: complex Morlet wavelet with bandwidth=1.5, center_frequency=1.0
# Morlet is standard for biological oscillation detection; good time-frequency tradeoff
wavelet_name = 'cmor1.5-1.0'
center_freq = pywt.central_frequency(wavelet_name)

# Test periods from 6h to 48h
# 6h minimum: twice the Nyquist-limited resolvable period
# 48h maximum: half the total duration; longer periods unreliable
periods_to_test = np.arange(6, 49, 0.5)
scales = center_freq * periods_to_test / sampling_interval

# CWT returns complex coefficients
coefficients, frequencies = pywt.cwt(signal_noisy, scales, wavelet_name, sampling_period=sampling_interval)

# Power: |coefficients|^2; represents energy at each time-frequency point
power_matrix = np.abs(coefficients) ** 2

# --- Extract dominant period over time ---
dominant_period_idx = np.argmax(power_matrix, axis=0)
dominant_period_over_time = periods_to_test[dominant_period_idx]

# --- Identify time windows with strong periodicity ---
# Threshold: mean + 2*SD of power; above-background significance
# This is a simple threshold; cone of influence correction improves edge accuracy
power_threshold = np.mean(power_matrix) + 2 * np.std(power_matrix)
max_power_over_time = np.max(power_matrix, axis=0)
significant_mask = max_power_over_time > power_threshold

n_significant = np.sum(significant_mask)
print(f'Timepoints with significant periodicity: {n_significant}/{n_timepoints}')

# --- Period transitions ---
# Detect when dominant period shifts
period_diff = np.abs(np.diff(dominant_period_over_time))
# Threshold of 3h: detects jumps between distinct oscillation regimes
transition_points = np.where(period_diff > 3.0)[0]
print(f'Period transition points (indices): {transition_points}')
if len(transition_points) > 0:
    print(f'Transition times (hours): {times[transition_points]}')

# --- Global wavelet spectrum (time-averaged power) ---
global_spectrum = np.mean(power_matrix, axis=1)
peak_indices, _ = find_peaks(global_spectrum, prominence=0.1 * np.max(global_spectrum))
print(f'\nGlobal spectrum peaks at periods: {periods_to_test[peak_indices]} hours')

# --- Multi-gene wavelet analysis ---
n_genes = 50
gene_results = []
for g in range(n_genes):
    if g < 15:
        # Genes with transient 12h oscillation (first half only)
        gene_signal = np.zeros(n_timepoints)
        gene_signal[:n_timepoints // 2] = 1.5 * np.sin(2 * np.pi * times[:n_timepoints // 2] / 12.0)
        gene_signal += np.random.normal(0, 0.4, n_timepoints)
    elif g < 30:
        # Genes with persistent 24h oscillation
        gene_signal = 2.0 * np.sin(2 * np.pi * times / 24.0) + np.random.normal(0, 0.4, n_timepoints)
    else:
        # Non-periodic genes
        gene_signal = np.random.normal(0, 0.5, n_timepoints)

    coeffs, _ = pywt.cwt(gene_signal, scales, wavelet_name, sampling_period=sampling_interval)
    gene_power = np.abs(coeffs) ** 2
    max_power = np.max(gene_power)
    dom_period = periods_to_test[np.argmax(np.mean(gene_power, axis=1))]

    gene_results.append({'gene': f'gene_{g}', 'max_power': max_power,
                          'dominant_period': dom_period})

# Simple power threshold for periodic gene detection
power_cutoffs = [r['max_power'] for r in gene_results]
# 70th percentile: separates periodic from non-periodic genes
cutoff_70 = np.percentile(power_cutoffs, 70)
periodic_genes = [r for r in gene_results if r['max_power'] > cutoff_70]
print(f'\nPeriodic genes detected (power > 70th percentile): {len(periodic_genes)}/{n_genes}')

# --- Visualization ---
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

axes[0, 0].plot(times, signal_noisy, 'steelblue', alpha=0.7, linewidth=0.8, label='Noisy')
axes[0, 0].plot(times, signal, 'red', linewidth=1.5, alpha=0.5, label='True')
axes[0, 0].axvline(48, color='black', linestyle='--', alpha=0.5, label='Period transition')
axes[0, 0].set_xlabel('Time (hours)')
axes[0, 0].set_ylabel('Expression')
axes[0, 0].set_title('Signal with time-varying periodicity')
axes[0, 0].legend(fontsize=8)

im = axes[0, 1].pcolormesh(times, periods_to_test, power_matrix, shading='auto', cmap='viridis')
axes[0, 1].axhline(12, color='white', linestyle='--', alpha=0.7, linewidth=0.8)
axes[0, 1].axhline(24, color='white', linestyle='--', alpha=0.7, linewidth=0.8)
axes[0, 1].axvline(48, color='red', linestyle='--', alpha=0.7, linewidth=0.8)
axes[0, 1].set_xlabel('Time (hours)')
axes[0, 1].set_ylabel('Period (hours)')
axes[0, 1].set_title('Wavelet scalogram')
axes[0, 1].invert_yaxis()
plt.colorbar(im, ax=axes[0, 1], label='Power')

axes[1, 0].plot(periods_to_test, global_spectrum, 'steelblue', linewidth=2)
for pi in peak_indices:
    axes[1, 0].axvline(periods_to_test[pi], color='red', linestyle='--', alpha=0.7)
    axes[1, 0].annotate(f'{periods_to_test[pi]:.0f}h',
                         xy=(periods_to_test[pi], global_spectrum[pi]),
                         xytext=(5, 5), textcoords='offset points', fontsize=9)
axes[1, 0].set_xlabel('Period (hours)')
axes[1, 0].set_ylabel('Mean power')
axes[1, 0].set_title('Global wavelet spectrum')

axes[1, 1].plot(times, dominant_period_over_time, 'steelblue', linewidth=1.5)
axes[1, 1].fill_between(times, dominant_period_over_time, alpha=0.2, color='steelblue')
axes[1, 1].axhline(12, color='red', linestyle=':', alpha=0.5, label='12h')
axes[1, 1].axhline(24, color='orange', linestyle=':', alpha=0.5, label='24h')
axes[1, 1].axvline(48, color='black', linestyle='--', alpha=0.5, label='Transition')
axes[1, 1].set_xlabel('Time (hours)')
axes[1, 1].set_ylabel('Dominant period (hours)')
axes[1, 1].set_title('Dominant period over time')
axes[1, 1].legend(fontsize=8)

plt.tight_layout()
plt.savefig('wavelet_periodicity_results.png', dpi=150, bbox_inches='tight')
plt.show()
print('\nPlot saved to wavelet_periodicity_results.png')
