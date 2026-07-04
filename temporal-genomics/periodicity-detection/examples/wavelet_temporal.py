# Reference: numpy 2.2+, scipy 1.15+, PyWavelets 1.8+, matplotlib 3.8+ | Verify API if version differs
# Transient / time-varying periodicity with a complex Morlet CWT, a cone-of-influence mask
# (edge artifacts), and Torrence & Compo (1998) AR(1) red-noise chi-square significance.
import os
import tempfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pywt
from scipy.stats import chi2
from scipy.signal import find_peaks

rng = np.random.default_rng(42)

# --- Signal whose period SHIFTS mid-record (global methods would blur this) ---
n = 200
times = np.linspace(0, 96, n)
dt = times[1] - times[0]                              # even sampling required by CWT
signal = np.where(times < 48,
                  2.0 * np.sin(2 * np.pi * times / 12.0),    # 0-48h: 12 h ultradian
                  1.5 * np.sin(2 * np.pi * times / 24.0))    # 48-96h: 24 h circadian
signal = signal + rng.normal(0, 0.4, n)

# --- Continuous Wavelet Transform (complex Morlet -> amplitude and phase) ---
wavelet = 'cmor1.5-1.0'
C = pywt.central_frequency(wavelet)                  # 1.0 for cmor1.5-1.0
periods_to_test = np.arange(6, 49, 0.5)
scales = C * periods_to_test / dt                    # scale = C * period / dt
assert np.allclose(pywt.scale2frequency(wavelet, scales) / dt, 1 / periods_to_test)
coeffs, _ = pywt.cwt(signal, scales, wavelet, sampling_period=dt)   # pass sampling_period
power = np.abs(coeffs) ** 2                           # power = |coefficients|^2

# --- Cone of influence: mask edge artifacts BEFORE reading ridges ---
# Morlet e-folding half-width ~ sqrt(2)*scale in time; with period = scale*dt (C=1), a period
# is inside the COI (unreliable) where sqrt(2)*period > distance-to-nearest-edge.
edge_dist = np.minimum(np.arange(n), np.arange(n)[::-1]) * dt
coi_max_period = edge_dist / np.sqrt(2)              # longest reliable period at each time
inside_coi = periods_to_test[:, None] > coi_max_period[None, :]
power_masked = np.where(inside_coi, np.nan, power)

# --- Torrence-Compo AR(1) red-noise significance (NOT mean + 2*SD) ---
alpha = float(np.corrcoef(signal[:-1], signal[1:])[0, 1])   # lag-1 autocorrelation
variance = signal.var(ddof=1)
# Theoretical red-noise background per period; local power/background ~ chi2 with 2 dof.
Pk = (1 - alpha**2) / (1 - 2 * alpha * np.cos(2 * np.pi * dt / periods_to_test) + alpha**2)
# pywt's Morlet power is not in variance units, so calibrate a single wavelet constant k
# from the ROBUST (median) noise level of the scalogram itself -- the median is dominated
# by the many noise cells, giving the expected white-noise power per unit variance.
k = np.nanmedian(power_masked / (variance * Pk[:, None]))
sig95 = variance * Pk * k * chi2.ppf(0.95, df=2) / 2   # 95% level, 2 dof (complex wavelet)
significant = power_masked > sig95[:, None]            # inside COI is NaN -> False
print(f'AR(1) lag-1 alpha={alpha:.3f}; significant time-period cells (COI-masked): {int(np.nansum(significant))}')

# --- Ridge: dominant period over time, COI-masked and above the red-noise contour ---
ridge = np.full(n, np.nan)
for j in range(n):
    col = np.where(significant[:, j], power_masked[:, j], np.nan)
    if np.any(np.isfinite(col)):
        ridge[j] = periods_to_test[np.nanargmax(col)]
early = np.nanmedian(ridge[times < 40])
late = np.nanmedian(ridge[times > 56])
print(f'Ridge period early (t<40h): {early:.1f}h | late (t>56h): {late:.1f}h  (expect ~12h -> ~24h)')

# --- Global wavelet spectrum (time-average, the wavelet analogue of an LS spectrum) ---
# Longest periods can be fully inside the COI (all-NaN row); average only finite cells.
finite_per_row = np.isfinite(power_masked).sum(axis=1)
global_spec = np.where(finite_per_row > 0, np.nansum(np.nan_to_num(power_masked), axis=1) / np.maximum(finite_per_row, 1), np.nan)
peaks, _ = find_peaks(global_spec, prominence=0.1 * np.nanmax(global_spec))
print(f'Global spectrum peaks at periods: {np.round(periods_to_test[peaks], 1)} h')

# --- Visualization (temp dir, removed at the end: no strays) ---
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
im = ax[0].pcolormesh(times, periods_to_test, power_masked, shading='auto', cmap='viridis')
ax[0].plot(times, np.clip(coi_max_period, periods_to_test.min(), periods_to_test.max()),
           color='white', ls='--', lw=1, label='COI')
ax[0].scatter(times, ridge, s=6, c='red', label='ridge')
ax[0].set_xlabel('Time (h)'); ax[0].set_ylabel('Period (h)'); ax[0].set_title('Scalogram (COI-masked)')
ax[0].invert_yaxis(); ax[0].legend(fontsize=8); plt.colorbar(im, ax=ax[0], label='Power')
ax[1].plot(global_spec, periods_to_test, 'steelblue', lw=2)
ax[1].scatter(global_spec[peaks], periods_to_test[peaks], c='red', zorder=3)
ax[1].set_xlabel('Mean power'); ax[1].set_ylabel('Period (h)'); ax[1].set_title('Global wavelet spectrum')
ax[1].invert_yaxis()
plt.tight_layout()
tmp = tempfile.mkdtemp()
out = os.path.join(tmp, 'wavelet_results.png')
plt.savefig(out, dpi=120); plt.close(fig)
os.remove(out); os.rmdir(tmp)
print('Figure rendered and cleaned up (temp dir removed)')
