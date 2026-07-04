---
name: bio-temporal-genomics-periodicity-detection
description: Discovers a periodic signal of UNKNOWN period in time-series omics data and puts a defensible significance on it, especially when sampling is IRREGULAR (dropped timepoints, pooled harvests) so FFT/Welch/JTK are invalid. Estimates the dominant period with Lomb-Scargle / generalized Lomb-Scargle (scipy, astropy), corroborates with autocorrelation, resolves transient/time-varying periodicity with the wavelet CWT (pywt), and screens genome-wide with false-alarm probabilities under BH FDR. Use when finding an oscillation whose period is not known a priori, analyzing cell-cycle or ultradian rhythms, or handling unevenly sampled time courses. Not for testing a KNOWN 24-hour rhythm (see temporal-genomics/circadian-rhythms).
tool_type: python
primary_tool: scipy
---

## Version Compatibility

Reference examples tested with: numpy 2.2+, scipy 1.15+ (lombscargle floating_mean present), astropy 8.0+, PyWavelets 1.8+, statsmodels 0.14+, matplotlib 3.8+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- scipy: `scipy.signal.lombscargle` changed in 1.17 (`precenter` deprecated -> removed 1.19; use `floating_mean=True` or pre-center). Check `scipy.__version__`.

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Periodicity Detection

**"Find a periodic pattern of unknown period in my time-series data"** -> Estimate the dominant period, attach a false-alarm probability, and (for transient signals) localize it in time. This is the complement of temporal-genomics/circadian-rhythms, which TESTS a known 24 h period; here the period is unknown and sampling is often uneven.

## Governing Principle

Unknown-period discovery is an ESTIMATION problem stacked on a DETECTION problem, and uneven sampling corrupts both. (1) Estimation: the period is a continuous quantity with a confidence region, not a yes/no on 24 h. (2) Detection: a periodogram peak is a random variable even under pure noise, and the MAX over a frequency grid is extreme-value-distributed, so significance is the False Alarm Probability (FAP) of that max, NOT the raw peak height. (3) Sampling design dominates: regular-grid methods (FFT, Welch, JTK, RAIN) assume equal spacing; genomics rarely delivers it (a timepoint fails QC, harvests are pooled, sampling is denser early). For uneven sampling the Lomb-Scargle least-squares periodogram replaces the FFT.

The single most important consequence: do NOT interpolate-then-FFT. Interpolation is a low-pass filter that injects spurious low-frequency power, suppresses power near the average Nyquist, and biases the whole spectrum red. Analyze the uneven series directly with Lomb-Scargle; interpolate only when a wavelet transform forces a grid, and flag it as an assumption.

## Method Selection

| Method | Sampling | Estimates | Significance model | Fails when |
|--------|----------|-----------|--------------------|-----------|
| Generalized Lomb-Scargle (astropy default, or scipy `floating_mean=True`) | Uneven OK | Global dominant period(s), amplitude, phase | Baluev analytic FAP (screen); bootstrap (few hits) | <2 cycles; strong harmonics; red noise mistaken for a peak |
| Classic Lomb-Scargle (scipy, pre-centered) | Uneven OK | Global dominant period | Permutation / Baluev via astropy | Nonzero-mean data if not centered; sparse points |
| Autocorrelation (statsmodels) | Even only | Fundamental period (coarse) | Bartlett bands (weak) | Uneven sampling; trend; poor resolution -> use as CHECK only |
| Wavelet CWT (pywt Morlet) | Even (interp small gaps) | Time-varying / transient period | Torrence-Compo AR(1) chi-square + COI mask | Edge claims (COI); needs a grid -> interpolation caveat |
| Welch PSD (scipy) | Even required | Smoothed global PSD | Segment-averaging variance reduction | Any gaps; long/close periods vs `nperseg` |
| Fisher's g-test (Wichert 2004) | Even only | Single dominant frequency, exact p | Exact analytic (white null) | Uneven sampling; red noise (white-null limitation) |

Decision spine: uneven times -> GLS + Baluev FAP -> BH FDR for a genome-wide screen; suspected transient -> CWT with COI + red-noise chi-square; even and want a quick classical p -> Fisher's g; ACF/Welch are even-sampling sanity checks, never the primary estimate. Re-verify the astropy FAP methods and scipy `floating_mean`/`normalize` semantics against the installed versions before relying on them.

## Lomb-Scargle Periodogram

**Goal:** Estimate the dominant period of an unevenly sampled series and attach a false-alarm probability.

**Approach:** Fit a sinusoid (plus a floating offset) at each trial frequency over a grid whose min/max and density are set from the record length, then read the peak period and its FAP. Prefer the generalized LS (floating mean) so a poorly-determined offset cannot leak into the amplitude.

Lomb-Scargle is NOT an FFT with jitter tolerance; it is a least-squares sinusoid fit at each frequency, invariant to time-origin, and (for even sampling) equivalent to the classical periodogram. Classical LS assumes zero-mean data; expression is always positive, so a nonzero mean leaks into the sinusoid and inflates power. Two correct fixes: pre-center (`y - y.mean()`, plug-in mean, old-school) or fit a per-frequency floating mean (generalized LS, the modern default).

### scipy: classic (pre-centered) vs generalized

```python
import numpy as np
from scipy.signal import lombscargle

# scipy uses ANGULAR frequency (rad/time): omega = 2*pi/period. Mixing ordinary
# frequency (1/period) is the #1 units bug -> periods off by 2*pi.
periods = np.linspace(6.0, 72.0, 2000)   # biologically plausible periods
angular = 2 * np.pi / periods

# Classic LS assumes zero mean: pre-center or power is inflated/distorted.
power_classic = lombscargle(times, values - values.mean(), angular, normalize=True)

# Generalized LS: fit a floating offset per frequency (Zechmeister & Kuerster).
# Preferred with sparse/uneven data; pass RAW values. floating_mean added in scipy
# 1.15; before that, pre-center. `precenter` is deprecated in 1.17, removed in 1.19.
power_gls = lombscargle(times, values, angular, normalize=True, floating_mean=True)

dominant_period = periods[np.argmax(power_gls)]
```

### astropy: generalized LS by default, with real FAP

**Goal:** Get a peak period plus a Baluev false-alarm probability, using a grid astropy sizes automatically.

**Approach:** astropy's LombScargle fits a floating mean by default; pass RAW values, let `autopower` build a correctly oversampled grid, then convert the peak power to a FAP.

```python
from astropy.timeseries import LombScargle

# fit_mean=True and center_data=True by DEFAULT -> generalized LS out of the box.
# Pass RAW values; pre-centering here is redundant. Uses ORDINARY frequency
# (cycles/time), NOT angular -> do not reuse a scipy omega-grid.
ls = LombScargle(times, values)          # add dy=sigma for heteroskedastic weighting
freq, power = ls.autopower(
    minimum_frequency=1 / 72.0,          # <= 1/record-span; need >=2 cycles to trust
    maximum_frequency=1 / 6.0,           # pseudo-Nyquist; uneven times can exceed 1/(2*dt_mean)
    samples_per_peak=10)                 # oversample: a peak has finite width ~1/T_span

best_period = 1 / freq[np.argmax(power)]

# FAP = P(noise produces a peak this high anywhere on the grid). It is GRID-DEPENDENT
# (rises with wider range / more oversampling), so always report the grid.
fap = ls.false_alarm_probability(power.max(), method='baluev')   # analytic, fast, screen-safe
levels = ls.false_alarm_level([0.1, 0.05, 0.01], method='baluev')
# method='bootstrap' (method_kwds={'n_bootstraps': 1000}) is most faithful but ~1000x cost
```

### Frequency-grid design (where peaks get missed)

The grid is not cosmetic: a bad grid silently loses signals.
- Minimum frequency comes from the record length: a period longer than ~the span cannot be claimed. Require >=2 full cycles (`f_min ~ 1/(T_span/2)`); 1 cycle cannot be told from a trend.
- Maximum frequency is a pseudo-Nyquist. There is no single Nyquist for irregular times; the average `1/(2*dt_mean)` is only a heuristic, and genuine unevenness permits probing far above it. Set `nyquist_factor` to a few, understanding it is a convention.
- Oversampling: a peak has finite width ~`1/T_span`; a grid coarser than a fraction of that steps over the true peak and under-reports its height and location. Use `samples_per_peak` 5-10 (astropy's `autopower` does this) rather than an ad-hoc `linspace(...,1000)`.
- Leakage / aliasing: gapped sampling has a window function that convolves the true spectrum, creating sidelobes and aliases at `1/(1/P +/- 1/P_sampling)`. Eyeball the window-function periodogram (LS of a constant signal at the same times) to see where sampling itself manufactures peaks.

## Autocorrelation (a sanity check, not an estimator)

**Goal:** Corroborate an LS-estimated period, not measure it.

**Approach:** Detrend, then look for a harmonic comb of ACF peaks at multiples of the period; treat a single band-crossing as weak evidence.

```python
from statsmodels.tsa.stattools import acf

# statsmodels assumes EVEN sampling and takes NO time vector -> invalid on the
# uneven data that motivates this skill. Detrend first: a trend gives a slowly
# decaying ACF that mimics periodicity.
vals_dt = values_even - np.polyval(np.polyfit(np.arange(len(values_even)), values_even, 1), np.arange(len(values_even)))
acf_vals, confint = acf(vals_dt, nlags=len(vals_dt) // 2, alpha=0.05)   # index 0 is lag-0 = 1.0
```

ACF is a WEAK period estimator: its resolution is quantized to the sampling step (period only to +/- one step), it requires even sampling, and trends fool it. A periodic signal repeats at lags `P, 2P, 3P`; the partial ACF (PACF) helps isolate the fundamental from harmonic echoes. Use it as a robustness box after detrending, never as the primary estimate.

## Wavelet CWT (transient / time-varying periodicity)

Global LS/Welch report one spectrum for the whole record and blur a signal that oscillates only early (cell-cycle synchrony decaying as cells desynchronize) or shifts period after a stimulus. The CWT resolves power in the time x period plane, at the cost of lower frequency precision and edge artifacts.

**Goal:** Map how the dominant period changes over time and mark which of that map is trustworthy and significant.

**Approach:** Transform with a complex Morlet, mask the cone of influence (edge artifacts), and test power against a Torrence-Compo AR(1) red-noise background rather than an ad-hoc `mean + 2*SD`.

```python
import pywt
import numpy as np
from scipy.stats import chi2

# Complex Morlet 'cmorB-C' (bandwidth B, center freq C): complex -> amplitude AND phase.
wavelet = 'cmor1.5-1.0'
C = pywt.central_frequency(wavelet)        # 1.0 for cmor1.5-1.0
dt = times[1] - times[0]                    # requires even sampling
periods_to_test = np.arange(6, 49, 0.5)
scales = C * periods_to_test / dt           # scale = C * period / dt ; verify below
# sanity: pywt.scale2frequency(wavelet, scales) / dt  ->  1/periods_to_test

# Pass sampling_period, else returned freqs are in per-sample units (period axis off by dt).
coeffs, freqs = pywt.cwt(signal, scales, wavelet, sampling_period=dt)
power = np.abs(coeffs) ** 2                  # power = |coefficients|^2
```

### Cone of influence (COI) - mask edge artifacts BEFORE reading ridges

Near each end of a finite record a wavelet overlaps the edge and its coefficients are computed against padding - edge artifacts, not signal. The COI widens for LONGER periods (larger scales reach farther from the edge). For the Morlet the e-folding half-width is ~`sqrt(2)*scale` in time; since `period = scale*dt/C` and C=1, a period is inside the COI (unreliable) where `sqrt(2)*period > distance-to-nearest-edge`. Grey out / mask the two triangular corners before ridge extraction or peak reading. Claiming a long-period oscillation that lives only in the first/last fraction of the record, inside the COI, is a classic false positive.

```python
n = len(signal)
edge_dist = np.minimum(np.arange(n), np.arange(n)[::-1]) * dt   # time to nearest edge
coi_max_period = edge_dist / np.sqrt(2)                          # longest reliable period per time
inside_coi = periods_to_test[:, None] > coi_max_period[None, :]
power_masked = np.where(inside_coi, np.nan, power)               # ignore masked cells downstream
```

### Red-noise significance (Torrence & Compo 1998), NOT mean + 2*SD

`mean(power) + 2*SD(power)` is indefensible: wavelet power under noise is not Gaussian and its variance changes with scale. Model the null as red noise - an AR(1) process with lag-1 autocorrelation `alpha` estimated from the data (white noise, alpha=0, is too permissive for intrinsically red / 1/f omics). The theoretical background is `P_k = (1-alpha^2)/(1 - 2*alpha*cos(2*pi*dt/period) + alpha^2)`; local power normalized by this background is chi-square with 2 dof (complex wavelet), so the 95% contour is `(noise level) * P_k * chi2(0.95, 2)/2`. Peaks poking through it are significant. Report the `alpha` assumed. pywt's Morlet power is not in variance units, so calibrate one wavelet constant `k` from the ROBUST (median) noise level of the scalogram - the median is dominated by noise cells, so it estimates the white-noise power per unit variance without the strong signal cells inflating it.

```python
alpha = np.corrcoef(signal[:-1], signal[1:])[0, 1]   # lag-1 autocorrelation
variance = signal.var(ddof=1)
Pk = (1 - alpha**2) / (1 - 2*alpha*np.cos(2*np.pi*dt/periods_to_test) + alpha**2)
k = np.nanmedian(power_masked / (variance * Pk[:, None]))   # robust wavelet power constant
sig95 = variance * Pk * k * chi2.ppf(0.95, df=2) / 2 # per-period 95% level, 2 dof
significant = power_masked > sig95[:, None]          # inside COI is NaN -> False
```

### Ridge extraction

Track `argmax` over scale at each time for the instantaneous dominant period, but only AFTER masking the COI and only for points above the red-noise contour; enforce continuity (a real ridge does not teleport between distant periods sample-to-sample). The global wavelet spectrum (time-average of power) is the wavelet analogue of a Fourier/LS spectrum, with its own chi-square test at reduced dof.

## Welch PSD (evenly sampled only)

```python
from scipy.signal import welch

# nperseg is the bias/variance knob: longer segments -> finer frequency resolution but
# fewer segments -> noisier PSD; shorter -> smoother but cannot resolve close/long periods.
# n//2 is a middling, defensible-but-arbitrary compromise; you cannot see a period longer
# than one segment. Welch CANNOT handle gaps -> feeding it interpolated data reintroduces
# interpolation bias. Even sampling only; otherwise use Lomb-Scargle.
freqs_w, psd = welch(values_even, fs=1/dt, nperseg=len(values_even)//2, detrend='constant')
periods_w = 1 / freqs_w[1:]
```

## Genome-Wide Screening

**Goal:** Turn per-gene FAPs into a genome-wide error rate without inflating the hit list with harmonics or trends.

**Approach:** Compute one FAP per gene, control FDR across genes, and guard against harmonic contamination and non-white noise.

```python
from statsmodels.stats.multitest import multipletests

# FAP is per-gene: thresholding 15,000 genes at FAP<0.01 yields ~150 false positives.
# Convert per-gene FAPs to q-values with Benjamini-Hochberg (assumes independence /
# positive dependence; gene-gene correlation makes it mildly conservative).
reject, qvals, _, _ = multipletests(fap_per_gene, method='fdr_bh')
n_periodic = int((qvals < 0.05).sum())
```

Screening subtleties:
- Harmonic contamination: a non-sinusoidal 24 h oscillation has real power at 12 h, 8 h, 6 h. The 12 h peak is a genuine harmonic, NOT an independent 12 h rhythm. A naive screen reports a phantom cohort of 12 h genes. Guard: check whether a putative P/2 peak co-occurs with a stronger peak at P; fit the fundamental plus harmonics jointly with astropy `LombScargle(t, y, nterms=k)` and attribute power correctly; distrust exact 2:1 period ratios.
- Permutation null validity: shuffling values across the fixed times destroys ALL temporal structure -> the null becomes WHITE noise. That is correct only if the alternative is "periodicity vs i.i.d. noise." Real omics noise is autocorrelated / red (1/f): a gene with a smooth trend or slow drift beats a white null and is falsely flagged periodic. On un-detrended, autocorrelated data a shuffle test is ANTI-CONSERVATIVE. Fixes: detrend first; use an AR(1) surrogate / block-bootstrap null that preserves short-range autocorrelation; or use the analytic Baluev FAP. State which null was used and what alternative it implies.
- Fisher's g-test (Wichert, Fokianos & Strimmer 2004): for EVENLY sampled series, g = (max periodogram ordinate)/(sum of ordinates) has a known exact distribution under the white-noise null, giving an exact per-gene p with no simulation - the canonical microarray cell-cycle screen. Even sampling only; shares the white-null limitation.
- Interpolate-then-FFT biases the spectrum red; never interpolate uneven data to run an even-sampling screen - use LS/GLS directly.
- Genuine oscillation vs 1/f vs trend: a real oscillation is a NARROW peak riding above the smooth 1/f background and recurring at harmonics; 1/f humps are broad and non-harmonic; a trend is indistinguishable from a very long period over <2 cycles. Require >=2-3 cycles, detrend, and test against a colored (AR(1)) null, not white.

## Common Errors

| Trap | Why it is wrong | Fix |
|------|-----------------|-----|
| Raw (nonzero-mean) values into `scipy.signal.lombscargle` | Classic LS assumes zero mean; the offset leaks into the sinusoid -> inflated power | `floating_mean=True` (GLS) OR pass `y - y.mean()`; `precenter` deprecated in 1.17 |
| Ordinary frequency (1/period) in scipy | scipy takes ANGULAR omega=2*pi/period; period off by 2*pi | `angular = 2*np.pi/periods`; use astropy for ordinary-frequency grids |
| Pre-centering then handing to astropy | astropy already fits the mean (`fit_mean=True`) -> redundant / conceptual muddle | Pass RAW values to astropy; it is GLS by default |
| Ad-hoc `np.linspace(f_min,f_max,1000)` grid | Too-coarse spacing steps over the finite-width peak | `autopower(samples_per_peak=5..10)` or space finer than `~1/T_span` |
| Reporting FAP without the grid | FAP grows with grid width / oversampling; numbers not comparable | Always report `[f_min,f_max]` and `samples_per_peak` |
| Interpolate-then-FFT/Welch on uneven data | Interpolation is a low-pass filter -> spurious red power, killed high-freq signal | Analyze the uneven series directly with LS/GLS |
| Wavelet scalogram with no COI mask | Edge coefficients are artifacts against padding; worse at long periods | Compute the COI (Morlet half-width `sqrt(2)*scale`), grey it out |
| `mean + 2*SD` wavelet significance | Wavelet power is not Gaussian; variance varies with scale | Torrence-Compo AR(1) background x chi2(0.95,2)/2 contour; report alpha |
| Shuffle-time permutation on un-detrended data | White null; a trend / red-noise gene beats it -> false "periodic" | Detrend first, or AR(1)/block-bootstrap surrogate, or Baluev FAP |
| ACF as the period estimator | Poor resolution; assumes even sampling; trend-fooled | Use ACF only as a corroborating check after detrending |
| Treating the 12 h (P/2) peak as an independent rhythm | It is the harmonic of a non-sinusoidal 24 h signal | Check co-occurrence with a stronger peak at P; fit `nterms>1` |
| Claiming a period longer than the record | <2 cycles cannot separate oscillation from trend | Require >=2 (ideally 3) cycles; cap `max_period <= T_span/2` |
| Assuming 24 h for cell cycle | Cell-cycle period is cell-type / condition dependent, usually != 24 h | Estimate the period; use circadian-rhythms only when 24 h is the hypothesis |

## References

- Lomb 1976. Least-squares frequency analysis of unequally spaced data. Astrophys Space Sci 39(2):447-462.
- Scargle 1982. Studies in astronomical time series analysis. II. Statistical aspects of spectral analysis of unevenly spaced data. Astrophys J 263:835-853.
- VanderPlas 2018. Understanding the Lomb-Scargle Periodogram. Astrophys J Suppl Ser 236(1):16.
- Baluev 2008. Assessing the statistical significance of periodogram peaks. Mon Not R Astron Soc 385(3):1279-1285.
- Torrence & Compo 1998. A Practical Guide to Wavelet Analysis. Bull Amer Meteor Soc 79(1):61-78.
- Wichert, Fokianos & Strimmer 2004. Identifying periodically expressed transcripts in microarray time series data. Bioinformatics 20(1):5-20.

## Related Skills

temporal-genomics/circadian-rhythms - Known-period (24 h) rhythm testing with cosinor and JTK_CYCLE
temporal-genomics/temporal-clustering - Group genes by periodicity characteristics
temporal-genomics/trajectory-modeling - Non-periodic trajectory fitting with GAMs
