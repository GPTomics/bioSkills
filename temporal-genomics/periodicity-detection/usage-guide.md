# Periodicity Detection - Usage Guide

## Overview

Discovers a periodic signal of UNKNOWN period in time-series omics data and attaches a defensible significance to it. This is the complement of temporal-genomics/circadian-rhythms, which tests a KNOWN 24 h period: here the period is unknown and must be estimated, and sampling is often irregular (dropped timepoints, pooled harvests, denser sampling early). Unknown-period discovery is an estimation problem (what period best explains the data, with a confidence region) stacked on a detection problem (is any peak larger than the noise floor of the MAX over a frequency grid). Uneven sampling corrupts both, which is why the generalized Lomb-Scargle periodogram, not the FFT, is the native tool. Autocorrelation is a corroborating sanity check; the wavelet CWT resolves transient or time-varying periodicity that global methods blur; Welch and Fisher's g apply only to evenly sampled data.

## Prerequisites

### Python
```bash
pip install numpy scipy astropy PyWavelets statsmodels matplotlib
```

### Data Requirements
- Time-series expression values with corresponding observation times
- At least 2 complete cycles of the shortest period of interest (1 cycle cannot be told from a trend)
- Lomb-Scargle: observation times may be irregularly spaced (its reason to exist)
- Wavelet / Welch / autocorrelation / Fisher's g: evenly spaced data (interpolate only small, roughly uniform gaps, and flag it as an assumption)
- scipy 1.17 changed `scipy.signal.lombscargle` (`precenter` deprecated, use `floating_mean=True`); check `scipy.__version__`

## Quick Start

Tell the AI agent what to search for:
- "Find periodic expression patterns in my unevenly sampled time-series data and estimate the period"
- "What is the dominant oscillation period in my gene expression time course, and is it significant?"
- "Detect transient periodicity that starts strong then fades in my developmental time series"
- "Screen my RNA-seq time course for periodic genes with FDR control"

## Example Prompts

### Unknown Period Discovery
> "I have 96 hours of gene expression sampled at irregular intervals. Estimate which genes are periodic and what their periods are, and give me a false-alarm probability for each."

> "Compute a generalized Lomb-Scargle periodogram for this gene and report the dominant period with a Baluev FAP, using an oversampled frequency grid."

### Cell Cycle and Ultradian Oscillations
> "Search for cell-cycle periodicity in my expression data. I expect a period near 18-24 h but do not know it exactly; do not assume 24 h."

> "Which genes oscillate with a period different from 24 h, and could any 12 h hits just be harmonics of a 24 h signal?"

### Transient Periodicity
> "Some genes oscillate early then stop as the population desynchronizes. Use a wavelet CWT, mask the cone of influence, and mark which time-period regions are significant against a red-noise background."

> "Apply a continuous wavelet transform and show me a scalogram of how the dominant period changes over time."

### Genome-Wide Screening
> "Screen all 15,000 expressed genes for significant periodicity with Lomb-Scargle and BH FDR correction; tell me which null you used."

> "Run a permutation-based periodicity test but detrend first so drift does not create false positives."

## What the Agent Will Do

1. Load the time series and inspect the sampling structure (even vs uneven, span, gaps)
2. Detrend if a monotonic trend or drift is present (trends create low-frequency power that mimics a long period)
3. Estimate the dominant period with the generalized Lomb-Scargle periodogram, sizing the frequency grid (min from record length, max at a pseudo-Nyquist, oversampled so peaks are not stepped over)
4. Attach a false-alarm probability (Baluev analytic for screens, bootstrap for a few hits), always reporting the grid because FAP is grid-dependent
5. Corroborate with autocorrelation (after detrending) as a sanity check, not the primary estimate
6. For suspected transient periodicity, run a wavelet CWT, mask the cone of influence, and test power against a Torrence-Compo AR(1) red-noise chi-square background instead of an ad-hoc threshold
7. For genome-wide screens, convert per-gene FAPs to q-values with Benjamini-Hochberg, guard against harmonic contamination (the 12 h harmonic of a 24 h signal), and distinguish genuine oscillation from 1/f noise and trend
8. Export significant periodic genes with estimated periods, FAP/q-values, and the null/grid used

## Tips

- Generalized Lomb-Scargle (floating mean) is the default: it fits a per-frequency offset so a poorly-determined mean cannot leak into the amplitude. astropy fits the mean by default (pass RAW values); with scipy set `floating_mean=True` or pre-center.
- scipy uses ANGULAR frequency (omega = 2*pi/period); astropy uses ordinary frequency (1/period). Do not reuse one grid in the other.
- Never interpolate-then-FFT uneven data: interpolation is a low-pass filter that injects spurious low-frequency power and biases the spectrum red. Analyze the uneven series directly.
- FAP is grid-dependent: it rises as the frequency range widens or oversampling increases, so two "FAP = 0.01" values over different grids are not comparable. Always report `[f_min, f_max]` and `samples_per_peak`.
- Require at least 2 (ideally 3) full cycles within the record; a period near the record length is indistinguishable from a trend.
- Wavelet edges are artifacts: mask the cone of influence (widens for longer periods) before reading ridges or peaks, and test against a red-noise (AR(1)) background, not `mean + 2*SD`.
- Shuffling timepoints gives a WHITE null; on autocorrelated or trended data it is anti-conservative. Detrend first, or use an AR(1) surrogate / block-bootstrap, or the analytic Baluev FAP.
- A 12 h peak may be the harmonic of a non-sinusoidal 24 h signal, not an independent rhythm; fit the fundamental plus harmonics (`nterms>1`) and distrust exact 2:1 period ratios.
- For cell-cycle analysis the period depends on cell type and growth conditions; estimate it rather than assuming 24 h.
- Permutation FDR is expensive genome-wide (1000 permutations x 15,000 genes); prefer analytic Baluev FAP for the screen and pre-filter to variable genes.

## Related Skills

temporal-genomics/circadian-rhythms - Known-period (24 h) rhythm testing with cosinor and JTK_CYCLE
temporal-genomics/temporal-clustering - Group genes by periodicity characteristics
temporal-genomics/trajectory-modeling - Non-periodic trajectory fitting with GAMs
