# Reference: numpy 2.2+, scipy 1.15+, astropy 8.0+, statsmodels 0.14+, matplotlib 3.8+ | Verify API if version differs
# Unknown-period discovery on IRREGULARLY sampled data with generalized Lomb-Scargle,
# a defensible false-alarm probability (Baluev), and genome-wide BH FDR control.
import os
import tempfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import lombscargle
from astropy.timeseries import LombScargle
from statsmodels.stats.multitest import multipletests

rng = np.random.default_rng(42)

# --- Irregularly sampled periodic signal (mimics dropped/pooled timepoints) ---
n_obs = 40
times = np.sort(rng.uniform(0, 96, n_obs))          # 96 h span, uneven spacing
true_period = 18.0                                   # non-circadian (fast-dividing cell cycle)
amplitude, mesor = 2.0, 8.0                          # nonzero mean -> classic-LS centering matters
signal = mesor + amplitude * np.sin(2 * np.pi * times / true_period) + rng.normal(0, 0.5, n_obs)

# --- scipy generalized LS (floating mean); ANGULAR frequency grid ---
# min_period ~ span/2 demands >=2 cycles; max_freq at 6 h ~ pseudo-Nyquist for ~2.4 h median spacing.
periods = np.linspace(6.0, 72.0, 2000)               # 2000 pts: finer than peak width ~1/T_span
angular = 2 * np.pi / periods
# floating_mean=True (scipy >=1.15) fits a per-frequency offset -> pass RAW values.
power_gls = lombscargle(times, signal, angular, normalize=True, floating_mean=True)
detected_period = periods[np.argmax(power_gls)]
print(f'True period: {true_period:.1f}h | scipy-GLS detected: {detected_period:.1f}h')

# --- astropy: generalized LS by default + Baluev analytic FAP (ORDINARY frequency) ---
ls = LombScargle(times, signal)                      # fit_mean=True default -> pass RAW values
freq, power = ls.autopower(minimum_frequency=1/72.0, maximum_frequency=1/6.0, samples_per_peak=10)
best_period = 1 / freq[np.argmax(power)]
# FAP is grid-dependent; report the grid [1/72, 1/6] cyc/h, samples_per_peak=10.
fap = ls.false_alarm_probability(power.max(), method='baluev')
levels = ls.false_alarm_level([0.1, 0.05, 0.01], method='baluev')
print(f'astropy-GLS best period: {best_period:.1f}h | Baluev FAP: {fap:.2e}')

# --- Permutation null: DETREND first so the null is not beaten by drift ---
# Shuffling times yields a WHITE null; on autocorrelated/trended data it is anti-conservative.
# Removing a linear trend before shuffling keeps the test honest for this drift-free demo.
trend = np.polyval(np.polyfit(times, signal, 1), times)
resid = signal - trend
n_perm = 500
null_max = np.array([lombscargle(times, rng.permutation(resid), angular, normalize=True).max()
                     for _ in range(n_perm)])
obs = lombscargle(times, resid, angular, normalize=True).max()
perm_p = float(np.mean(null_max >= obs))
print(f'Permutation p (detrended, white null): {perm_p:.4f}')

# --- Genome-wide screen: per-gene Baluev FAP -> BH FDR across genes ---
n_genes, n_periodic = 200, 40
faps, dom_periods = [], []
for g in range(n_genes):
    if g < n_periodic:
        p = rng.uniform(10, 30)                       # diverse unknown periods
        y = rng.uniform(1, 3) * np.sin(2 * np.pi * times / p) + rng.normal(0, 0.5, n_obs)
    else:
        y = rng.normal(0, 0.5, n_obs)                 # aperiodic
    g_ls = LombScargle(times, y)
    g_freq, g_pow = g_ls.autopower(minimum_frequency=1/72.0, maximum_frequency=1/6.0, samples_per_peak=10)
    faps.append(g_ls.false_alarm_probability(g_pow.max(), method='baluev'))
    dom_periods.append(1 / g_freq[np.argmax(g_pow)])
faps = np.clip(np.array(faps), 1e-300, 1.0)           # Baluev bound can underflow to 0
reject, qvals, _, _ = multipletests(faps, method='fdr_bh')   # BH across genes
n_hit = int((qvals < 0.05).sum())
tp = int((qvals[:n_periodic] < 0.05).sum())
print(f'Genome-wide: {n_hit}/{n_genes} periodic at q<0.05 | true positives {tp}/{n_periodic}')

# --- Visualization (written to a temp dir, then removed: leaves no strays) ---
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
ax[0].plot(1 / freq, power, 'steelblue', lw=1.2)
ax[0].axvline(true_period, color='red', ls='--', alpha=0.6, label=f'true {true_period}h')
for lv, lab in zip(levels, ['FAP 0.1', 'FAP 0.05', 'FAP 0.01']):
    ax[0].axhline(lv, color='gray', ls=':', alpha=0.6, label=lab)
ax[0].set_xlabel('Period (h)'); ax[0].set_ylabel('Normalized power')
ax[0].set_title(f'astropy GLS (peak {best_period:.1f}h)'); ax[0].legend(fontsize=8); ax[0].invert_xaxis()
mask = qvals < 0.05
ax[1].scatter(np.array(dom_periods)[~mask], -np.log10(qvals[~mask]), c='gray', s=18, alpha=0.5, label='ns')
ax[1].scatter(np.array(dom_periods)[mask], -np.log10(qvals[mask]), c='coral', s=28, label='q<0.05')
ax[1].axhline(-np.log10(0.05), color='black', ls='--', alpha=0.5)
ax[1].set_xlabel('Dominant period (h)'); ax[1].set_ylabel('-log10(q)')
ax[1].set_title(f'Genome-wide screen ({n_hit} hits)'); ax[1].legend(fontsize=8)
plt.tight_layout()
tmp = tempfile.mkdtemp()
out = os.path.join(tmp, 'lomb_scargle_results.png')
plt.savefig(out, dpi=120); plt.close(fig)
os.remove(out); os.rmdir(tmp)
print(f'Figure rendered and cleaned up (temp dir removed)')
