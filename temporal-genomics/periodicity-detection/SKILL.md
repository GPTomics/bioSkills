---
name: bio-temporal-genomics-periodicity-detection
description: Discovers periodic signals of unknown period in time-series omics data using Lomb-Scargle periodograms (scipy), autocorrelation, and wavelet time-frequency decomposition (pywt). Identifies dominant frequencies, handles irregularly sampled data, and detects transient periodicity. Use when searching for periodic patterns of unknown period length, analyzing cell cycle oscillations, or processing unevenly spaced time-series. Not for testing known 24-hour rhythms (see temporal-genomics/circadian-rhythms).
tool_type: python
primary_tool: scipy
---

# Periodicity Detection

Discovers periodic signals of unknown frequency in time-series omics data. Handles irregular sampling, identifies dominant oscillation periods, and detects transient or time-varying periodicity through spectral and time-frequency methods.

## Core Workflow

1. Prepare time-series expression data (handle missing values, detrend if needed)
2. Compute frequency spectrum (Lomb-Scargle, Welch, or autocorrelation)
3. Identify significant spectral peaks
4. Assess statistical significance (false alarm probability, permutation FDR)
5. For transient periodicity, apply wavelet time-frequency decomposition

## Lomb-Scargle Periodogram (scipy)

Handles irregularly sampled time series without interpolation. Standard method for astronomical and biological time-series with uneven sampling.

### Basic Lomb-Scargle

```python
import numpy as np
from scipy.signal import lombscargle

# times: observation times (may be unevenly spaced)
# values: expression values at those times
# Subtract mean for proper normalization
values_centered = values - np.mean(values)

# Frequency grid: test periods from 4h to 72h
# min_period=4h: Nyquist limit for ~2h sampling; shorter periods unresolvable
# max_period=72h: longer than this needs more data (at least 2 full cycles)
min_freq = 2 * np.pi / 72.0
max_freq = 2 * np.pi / 4.0
# 1000 frequencies: fine grid for precise period estimation
freqs = np.linspace(min_freq, max_freq, 1000)

# Angular frequencies for scipy.signal.lombscargle
power = lombscargle(times, values_centered, freqs, normalize=True)

# Convert to periods for interpretation
periods = 2 * np.pi / freqs
```

### Peak Detection and Significance

```python
from scipy.signal import find_peaks

# Find spectral peaks
# prominence=0.1: minimum peak prominence above local background
peaks, properties = find_peaks(power, prominence=0.1)

# Dominant period
dominant_idx = np.argmax(power)
dominant_period = periods[dominant_idx]
dominant_power = power[dominant_idx]
```

### False Alarm Probability (FAP)

```python
from astropy.timeseries import LombScargle

# astropy provides built-in FAP calculation
ls = LombScargle(times, values_centered)
frequency, power_astropy = ls.autopower(
    minimum_frequency=1 / 72.0,
    maximum_frequency=1 / 4.0
)

# FAP: probability that noise alone produces a peak this high
# FAP < 0.01: strong evidence of true periodicity
fap = ls.false_alarm_probability(power_astropy.max())

# FAP at multiple significance levels
fap_levels = ls.false_alarm_level([0.01, 0.05, 0.10])
```

### Genome-Wide Periodicity Screening

```python
from statsmodels.stats.multitest import multipletests

# Screen all genes; collect FAP for each
fap_values = []
dominant_periods = []
for gene_idx in range(expr_mat.shape[0]):
    vals = expr_mat[gene_idx] - np.mean(expr_mat[gene_idx])
    power = lombscargle(times, vals, freqs, normalize=True)
    dominant_periods.append(periods[np.argmax(power)])
    # Approximate FAP using Baluev (2008) method via astropy
    ls = LombScargle(times, vals)
    _, pwr = ls.autopower(minimum_frequency=1 / 72.0, maximum_frequency=1 / 4.0)
    fap_values.append(ls.false_alarm_probability(pwr.max()))

# BH FDR correction for genome-wide multiple testing
reject, qvals, _, _ = multipletests(fap_values, method='fdr_bh')
# q < 0.05: standard FDR threshold
n_periodic = np.sum(qvals < 0.05)
```

## Autocorrelation

Detects periodicity by measuring self-similarity at increasing time lags.

```python
from statsmodels.tsa.stattools import acf

# nlags: test up to 2x expected maximum period (in sampling intervals)
# For 4h sampling testing up to 72h: nlags = 72/4 = 18
acf_values, confint = acf(values, nlags=18, alpha=0.05)

# Periodic signal produces peaks at multiples of the period
# Find first significant peak after lag 0
# Peaks above upper confidence bound indicate significant autocorrelation
peak_lags, _ = find_peaks(acf_values[1:])
significant_peaks = peak_lags[acf_values[peak_lags + 1] > confint[peak_lags + 1, 1]]
if len(significant_peaks) > 0:
    # First significant peak lag corresponds to the fundamental period
    fundamental_period = (significant_peaks[0] + 1) * sampling_interval
```

## Wavelet Time-Frequency Decomposition (pywt)

Detects transient or time-varying periodicity that global methods miss.

### Continuous Wavelet Transform (CWT)

```python
import pywt
import numpy as np

# Morlet wavelet: standard for biological oscillation detection
# Balances time and frequency resolution
wavelet = 'cmor1.5-1.0'

# Scales correspond to periods being tested
# scales = center_freq * period * sampling_rate
# For 4h sampling, test periods 8-72h
sampling_rate = 1 / 4.0
center_freq = pywt.central_frequency(wavelet)
periods_to_test = np.arange(8, 73, 1)
scales = center_freq * periods_to_test * sampling_rate

# CWT returns complex coefficients; power = |coefficients|^2
coefficients, frequencies = pywt.cwt(values, scales, wavelet, sampling_period=4.0)
power_matrix = np.abs(coefficients) ** 2
```

### Scalogram Visualization

```python
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(10, 6))
im = ax.pcolormesh(times, periods_to_test, power_matrix, shading='auto', cmap='viridis')
ax.set_xlabel('Time (hours)')
ax.set_ylabel('Period (hours)')
ax.set_title('Wavelet scalogram')
ax.invert_yaxis()
plt.colorbar(im, ax=ax, label='Power')
```

### Detect Time-Varying Periodicity

```python
# Ridge extraction: track dominant period over time
dominant_period_over_time = periods_to_test[np.argmax(power_matrix, axis=0)]

# Identify time windows where periodicity is strong
# Threshold: mean + 2*SD of power; captures peaks above background noise
power_threshold = np.mean(power_matrix) + 2 * np.std(power_matrix)
significant_mask = np.max(power_matrix, axis=0) > power_threshold
```

## Welch Periodogram (Evenly Sampled)

For regularly sampled data, Welch's method provides smooth power spectral density estimates.

```python
from scipy.signal import welch

# nperseg: segment length; controls frequency resolution vs variance tradeoff
# Longer segments = finer frequency resolution but more variance
# nperseg = n_timepoints // 2: standard compromise
fs = 1 / sampling_interval
freqs_welch, psd = welch(values, fs=fs, nperseg=len(values) // 2)

# Convert to periods
periods_welch = 1 / freqs_welch[1:]
psd_trimmed = psd[1:]
```

## Permutation-Based FDR

For genome-wide periodicity screening without parametric assumptions.

```python
# Permutation test: shuffle timepoints to generate null distribution
# n_perm=1000: standard for genome-wide testing; increase to 10000 for publication
n_perm = 1000
null_max_powers = np.zeros(n_perm)
for perm in range(n_perm):
    shuffled = np.random.permutation(values)
    null_power = lombscargle(times, shuffled - np.mean(shuffled), freqs, normalize=True)
    null_max_powers[perm] = np.max(null_power)

# Empirical p-value: fraction of null peaks >= observed peak
empirical_p = np.mean(null_max_powers >= dominant_power)
```

## Method Selection

| Method | Sampling | Detects | Limitations |
|--------|----------|---------|-------------|
| Lomb-Scargle | Uneven OK | Global periodicity | No time localization |
| Autocorrelation | Even preferred | Fundamental period | Low frequency resolution |
| Wavelet CWT | Even preferred | Transient periodicity | Lower frequency precision |
| Welch PSD | Even required | Global PSD | Cannot handle gaps |

## Parameter Guide

| Parameter | Typical Value | Rationale |
|-----------|---------------|-----------|
| Min period | 2x sampling interval | Nyquist criterion |
| Max period | Total duration / 2 | Need at least 2 full cycles |
| FAP threshold | < 0.01 | Conservative for individual genes |
| FDR threshold | q < 0.05 | Standard for genome-wide screening |
| Wavelet | cmor1.5-1.0 (Morlet) | Standard for biological oscillations |
| Permutations | 1000-10000 | 1000 for screening, 10000 for publication |

## Related Skills

- circadian-rhythms - Known-period rhythm testing with cosinor and JTK_CYCLE
- temporal-clustering - Group genes by periodicity characteristics
- trajectory-modeling - Non-periodic trajectory fitting with GAMs
