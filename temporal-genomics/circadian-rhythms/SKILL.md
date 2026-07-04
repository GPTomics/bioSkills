---
name: bio-temporal-genomics-circadian-rhythms
description: Tests and estimates rhythmicity at a PRE-SPECIFIED period (canonically 24h) in time-series omics using cosinor regression (CosinorPy), JTK_CYCLE/ARSER/Lomb-Scargle meta-analysis (MetaCycle meta2d), and non-parametric tests for asymmetric waveforms (RAIN, DiscoRhythm); estimates phase (acrophase), amplitude, and MESOR, and controls FDR with an effect-size (rAMP) filter against over-detection. Use when testing for 24-hour or other known-period oscillations in a single condition (circadian, feeding-fasting, or light-dark experiments) and estimating their phase/amplitude. Not for unknown-period discovery (see temporal-genomics/periodicity-detection) or comparing rhythms between conditions (see temporal-genomics/differential-rhythmicity).
tool_type: mixed
primary_tool: CosinorPy
---

## Version Compatibility

Reference examples tested with: CosinorPy 3.1 (requires numpy<2.0 - v3.1 calls the removed `np.round_`), pandas 2.2+, statsmodels 0.14+, MetaCycle 1.2+, RAIN 1.x (Bioconductor), DiscoRhythm 1.x (Bioconductor).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

CosinorPy 3.1 imports as `from CosinorPy import cosinor, cosinor1, file_parser` (capitalized package, lowercase submodules); older 0.x/1.x releases used the lowercase `cosinorpy` package name.

# Known-Period Rhythm Testing

## Governing principle: known-period TESTING, not unknown-period DISCOVERY

This skill answers "at THIS period (usually 24h), is a feature rhythmic, and what are its phase, amplitude, and MESOR?" - a hypothesis test plus parameter estimation at a period the analyst specifies. That is categorically different from asking "what period does this feature have?", which is a spectral DISCOVERY search (Lomb-Scargle periodogram, wavelets, FFT) handled by temporal-genomics/periodicity-detection. Conflating them is the field's most common conceptual error: opening a wide period window turns a test into a search and inflates false positives, because structured noise can always be "fit" better at SOME period in a broad window.

The load-bearing consequence: temporal conclusions are dominated by SAMPLING DESIGN, not by the algorithm. Nyquist (>=2 samples/cycle) is a mathematical FLOOR that only prevents aliasing - it gives zero robustness to noise and no ability to estimate phase/amplitude. Real detection needs >=6 (ideally 8-12) samples/cycle AND >=2 full cycles. Resolving that a period exists < estimating its phase < estimating its amplitude, in ascending sampling demand. A design good enough to say "yes, 24h" is usually too thin to trust its phase and far too thin to trust its amplitude.

## Core Workflow

1. Declare the light regime (LD=entrained, use ZT; DD=free-running, use CT) and the period window (entrained: fix `minper=maxper=24`; free-running: allow ~22-26h for tau != 24h)
2. Prepare the time-series matrix (features x timepoints), decide log-vs-linear scale and any detrending BEFORE testing
3. Fit cosinor models or apply rhythmicity tests at the specified period
4. Extract parameters: amplitude, relative amplitude (rAMP), phase (acrophase), MESOR, p-value
5. Control FDR (BH), then apply an EFFECT-SIZE filter (rAMP / fold-change) - significance alone over-detects
6. For between-condition questions, fit a differential-rhythmicity model (never intersect two separate rhythm lists)

## Method Selection (which to pick and why)

| Method | Pick when | Mechanism | Fails / caveat |
|--------|-----------|-----------|----------------|
| Cosinor (single-component) | Sinusoidal waveform; uneven/sparse/non-integer sampling; CIs on phase/amplitude/MESOR are needed; substrate for differential rhythmicity | OLS of expression on a fixed `cos`/`sin` basis at period T (linear regression); rhythmicity = zero-amplitude F-test | Miscalls asymmetric/spiky waveforms (a fast-rise/slow-decay pulse) as arrhythmic; needs a variance-stabilizing transform for count data |
| Cosinor (multi-component) | Visibly non-sinusoidal shape AND dense sampling AND a biological reason (e.g. a known 12h "12h-clock" transcript) | Adds 12h (`n_components=2`), 8h (`=3`) harmonics; joint zero-amplitude F-test | Each harmonic costs 2 df; with 6-8 pts/cycle a 3-component model is near-saturated and fits noise - use AIC/BIC or automatic model selection |
| JTK_CYCLE | Evenly sampled at integer-hour intervals; robust rank-based test; genome-scale speed | Correlates the series against reference cosines of all phases (Jonckheere-Terpstra + Kendall tau); best phase = matched reference | Requires EVEN integer sampling, no gaps; with few timepoints/1 replicate the tau null is DISCRETE so p-values are quantized and ANTI-conservative (source of "everything is rhythmic") |
| eJTK / BooteJTK | Short/sparse or few-replicate series where JTK p-values are untrustworthy; asymmetric/spiky waveforms | Empirical (permutation/Gamma) null restores calibration; asymmetric reference-waveform library; BooteJTK adds replicate bootstrap + variance shrinkage | Slower; still cannot fully fix temporal autocorrelation |
| ARSER | Non-sinusoidal short series; combines time- and frequency-domain info | Estimates period from an AUTOregressive spectrum, then harmonic regression | Requires EVEN sampling, no missing values, no replicate structure; AR order unstable on very short/noisy series; wants denser sampling than JTK |
| RAIN | ASYMMETRIC waveforms (fast induction / slow decay); distribution-free | Umbrella/Mack-Wolfe (Jonckheere-Terpstra) test with SEPARATE rising and falling limbs | LOWER power than cosinor/JTK for genuinely symmetric sinusoids; gives a coarse phase/peak-shape, not clean amplitude CIs |
| MetaCycle meta2d | A robust consensus RANK across methods is wanted on a standard even design | Runs a subset of {ARS,JTK,LS}, combines p by Fisher's method -> `meta2d_pvalue` (BH -> `meta2d_BH.Q`), averages period, circular-averages phase | Fisher assumes INDEPENDENT p; ARS/JTK/LS on the same data are correlated, so `meta2d_BH.Q` is NOT a literal FDR (read it as a rank aid); `analysisStrategy='auto'` SILENTLY drops ARS/JTK on uneven/replicated data (may run LS only); averaged period is meaningless when methods disagree |

For between-condition comparison, see temporal-genomics/differential-rhythmicity - none of the single-condition tests above answer it correctly.

## Design constraints (upstream of any method; non-negotiable)

- >=2 full cycles (48h circadian minimum; 3 cycles / 72h improves power and reveals damping). One cycle cannot distinguish an oscillation from a monotone trend or a single transient.
- >=6, ideally 8-12+, samples/cycle. 2-4h spacing is standard; 1-2h is needed to resolve waveform shape or fast harmonics. Sparse designs are exactly where JTK's calibration fails.
- >=2-3 biological replicates/timepoint. Single-replicate designs cripple FDR calibration (no within-timepoint variance; empirical-null/bootstrap corrections cannot work). Replication in TIME and AT a timepoint buy different things - do not trade all of one for the other.
- Harvest-ORDER confound (the silent killer): collecting/extracting/sequencing timepoints in temporal order aliases any drift (reagent lots, RIN, lane position) PERFECTLY onto ZT and manufactures spurious 24h rhythms. No rhythmicity test detects this. Fix by DESIGN: randomize processing order, balance replicates across batches, model batch as a covariate (trivial in a limma/DESeq2 design). It cannot be repaired analytically because batch and the rhythm are the same axis.

## CosinorPy (Python)

**Goal:** Test each feature for rhythmicity at a known period and estimate amplitude, relative amplitude, acrophase, and MESOR with FDR control.

**Approach:** Fit cosine curves per feature with `fit_group` (batch), use its built-in BH `q` column (or recompute BH over a chosen correction set), then filter on both q and relative amplitude.

### Single- and multi-component fit

Fits `y = M + A*cos(2*pi*t/T + phi)` where M = MESOR (rhythm-adjusted midline, NOT the arithmetic mean unless sampling is balanced), A = amplitude, phi = acrophase stored as `atan2(-gamma, beta)` (usually negative).

```python
from CosinorPy import cosinor, cosinor1, file_parser

df = file_parser.read_csv('expression_timecourse.csv')  # long format: columns x (time), y (value), test (feature id)

# Single-component (sinusoidal). period=24: standard circadian period in hours.
# fit_me returns a 5-tuple: (results, statistics, rhythm_params, X_test, Y_fit_test).
single = cosinor.fit_me(df[df['test'] == 'Arntl']['x'].values,
                        df[df['test'] == 'Arntl']['y'].values,
                        period=24, n_components=1)

# Multi-component adds harmonics for non-sinusoidal shape; add ONLY with dense sampling + a biological reason.
# fit_me takes a SINGLE n_components (an int); n_components=2 adds one 12h harmonic to the 24h fundamental.
two_comp = cosinor.fit_me(df[df['test'] == 'Dbp']['x'].values,
                          df[df['test'] == 'Dbp']['y'].values,
                          period=24, n_components=2)

# To let CosinorPy PICK the harmonic order by information criterion, fit a range with fit_group over a
# candidate list, then select per feature with get_best_models (do not pass a list to fit_me).
group_multi = cosinor.fit_group(df, period=24, n_components=[1, 2, 3], plot=False)
best_models = cosinor.get_best_models(df, group_multi, n_components=[1, 2, 3])
```

### Batch analysis with built-in q-values

**Goal:** Score every feature genome-wide and keep confident, high-amplitude oscillators.

**Approach:** `fit_group` returns per-feature statistics INCLUDING a BH-adjusted `q` column; add an rAMP effect-size filter on top of q.

```python
import numpy as np
from statsmodels.stats.multitest import multipletests

# fit_group returns columns: test, period, n_components, p, q, p_reject, q_reject, RSS, R2, R2_adj,
# log-likelihood, amplitude, acrophase, mesor, peaks, heights, troughs, heights2, ME, resid_SE.
results = cosinor.fit_group(df, period=24, n_components=1, plot=False)

# 'q' is already BH-adjusted across the fitted group. Recompute BH only if the correction SET should differ
# (e.g. exclude non-expressed features first). Default multipletests method is Holm-Sidak, so pass fdr_bh explicitly.
valid = results['p'].notna()
results.loc[valid, 'q_bh'] = multipletests(results.loc[valid, 'p'], method='fdr_bh')[1]

# rAMP = amplitude / MESOR normalizes out expression level so calls are comparable across features.
# rAMP > 0.1 (>=10% of baseline) is a conventional biological-relevance floor - sweep it, do not treat as law.
results['rAMP'] = results['amplitude'] / results['mesor']
rhythmic = results[(results['q'] < 0.05) & (results['rAMP'] > 0.1)]
```

### Population-mean cosinor (replicated / multi-subject)

**Goal:** Get group-level amplitude/phase with CIs that propagate BETWEEN-subject variance, instead of pseudoreplicating.

**Approach:** Fit one cosinor per subject and combine the estimates - pooling all subjects' points into one fit understates uncertainty.

```python
# cosinor1.population_fit_cosinor returns a DICT with keys: test, names, values, means, confint (nested amp/acr/mesor CIs),
# p_value, p_amp, p_acr, p_mesor (all underscore; e.g. pop['confint']['amp'], pop['p_amp']).
pop = cosinor1.population_fit_cosinor(subject_df, period=24, plot_on=False)
# cosinor1.population_fit_group(df, period=24) batches this across groups; cosinor1.population_test_cosinor_pairs
# compares two populations' rhythms (a differential-rhythmicity test on replicated data).
```

Convert acrophase to peak-hour with `peak_h = (-acrophase) * T / (2*pi) % T`; sanity-check against a known clock gene (mouse liver Arntl/Bmal1 peaks ~CT22-0, Nr1d1 ~CT4-6, Dbp ~CT8-10).

## MetaCycle meta2d (R)

**Goal:** Produce a robust consensus rhythmicity rank on an evenly sampled design.

**Approach:** Run meta2d over {JTK,ARS,LS}; read `meta2d_BH.Q` as a ranking aid (Fisher over correlated nulls, not a literal FDR), and distrust the averaged period/phase when constituents disagree.

```r
library(MetaCycle)
# minper=maxper=24 for entrained (LD) data; 22-26 for free-running (DD) where tau != 24h.
# timepoints must match column order. ARS/JTK need EVEN integer sampling with no missing values / no replicates;
# analysisStrategy='auto' silently drops ineligible methods (may leave LS only) - check the per-method columns.
# timepoints span 0-68h at 4h resolution: >=2 full 24h cycles at ~6 samples/cycle (the design floor this skill sets).
meta2d(infile = 'expression_matrix.csv', filestyle = 'csv', outdir = 'metaout',
       timepoints = seq(0, 68, by = 4), cycMethod = c('JTK', 'ARS', 'LS'),
       minper = 24, maxper = 24, outputFile = TRUE, outRawData = FALSE)

res <- read.csv('metaout/meta2d_expression_matrix.csv')
# meta2d_pvalue (Fisher-combined), meta2d_BH.Q, meta2d_period, meta2d_phase (hours from ZT0, peak time),
# meta2d_Base (baseline/MESOR), meta2d_AMP, meta2d_rAMP (= AMP/Base). Filter on rank AND relative amplitude.
rhythmic <- res[res$meta2d_BH.Q < 0.05 & res$meta2d_rAMP > 0.1, ]
```

## RAIN (R/Bioconductor)

**Goal:** Detect ASYMMETRIC waveforms (fast induction, slow decay) that cosinor/JTK miss.

**Approach:** Transpose to one-row-per-timepoint, declare replicate count, adjust p for multiple testing.

```r
library(rain)
# x needs ONE ROW PER TIMEPOINT (transpose a features x timepoints matrix). deltat = sampling interval (h).
# nr.series = replicates per timepoint (interleaved r1t1,r2t1,r1t2,...). method='independent' vs 'longitudinal'
# sets replicate handling. peak.border controls the allowed rising-fraction (asymmetry) window.
res <- rain(t(expression_mat), period = 24, deltat = 4, nr.series = 2, method = 'independent')
res$q <- p.adjust(res$pVal, method = 'BH')  # output columns: pVal, phase, peak.shape, period
rhythmic <- res[res$q < 0.05, ]
```

## DiscoRhythm (R/Bioconductor)

**Goal:** Run Cosinor/JTK/LS/ARS under one interface with built-in QC/PCA (scripted or Shiny).

```r
library(DiscoRhythm)
se <- discoGetSimu(TRUE)                                   # bundled demo SummarizedExperiment
disco <- discoBatch(se, osc_method = 'CS', report = NULL, osc_period = 24)  # osc_method='CS'=Cosinor; report=NULL skips the HTML report
```

Comparing rhythms BETWEEN conditions (differential rhythmicity: gain/loss/phase-shift/amplitude-change with LimoRhyde/dryR/compareRhythms, and the detect-then-Venn anti-pattern) is a distinct analysis - see temporal-genomics/differential-rhythmicity. Do NOT infer "genes that lost rhythm in the KO" by subtracting two independently thresholded single-condition rhythm lists.

## Common Errors (trap -> fix)

| Trap | Fix |
|------|-----|
| Opening a wide period window on a known-period test | Fix `minper=maxper=24` (entrained) or 22-26h (free-running); a wide window is discovery, not testing, and inflates false positives |
| Trusting JTK p-values / BH-Q from a single-replicate sparse design | Expect anti-conservative, quantized p-values; use eJTK/BooteJTK (empirical/bootstrap null) and inspect the genome-wide p-value HISTOGRAM before believing FDR |
| Reading `meta2d_BH.Q` as a literal FDR | Fisher integration over correlated ARS/JTK/LS p-values is not calibrated; use it as a consensus RANK and distrust averaged period/phase when methods disagree |
| Calling reduced BULK amplitude "arrhythmic" | Ensemble amplitude damps from cell DESYNCHRONY too; report "reduced ensemble amplitude" and use single-cell or imaging assays to separate loss-of-rhythm vs loss-of-synchrony |
| Claiming an "endogenous circadian rhythm" from LD data | LD rhythms can be light/feeding-DRIVEN (masking); endogeneity requires free-running (DD/constant) conditions. Diurnal != circadian. Use ZT for entrained, CT for free-running |
| Claiming a rhythm is "clock-CONTROLLED" from wild-type data alone | Persistence in DD proves endogeneity, not clock control; genetic dependence needs a clock-gene perturbation (compare WT vs clock-mutant, see temporal-genomics/differential-rhythmicity) |
| Mixing phase units/conventions (radians vs hours, +phi vs -phi, ZT vs CT) | State the convention; convert CosinorPy acrophase via `peak_h = (-acrophase)*T/(2*pi) % T`; sanity-check against a known clock gene's phase |
| Ranking features by RAW amplitude across the genome | Raw amplitude scales with expression and normalization; use relative amplitude (AMP/MESOR) or peak-to-trough fold-change for cross-feature comparison and the amplitude filter |
| Reporting significant rhythms with NO effect-size filter | Significance alone over-detects (Laloum 2020); add an rAMP/fold-change cutoff and report the amplitude DISTRIBUTION of the hit list, not just the count |
| Trusting phase/amplitude POINT estimates for near-threshold features | Estimation is unreliable where detection is marginal; interpret parameters only for confidently rhythmic features |
| Overfitting with `n_components=3` on 6-8 points/cycle | Harmonics cost 2 df each; use AIC/BIC or automatic model selection; add harmonics only with dense sampling and a biological reason |
| Harvest-order drift confounded with ZT | Randomize PROCESSING order, balance replicates across batches, model batch as a covariate; no rhythmicity test detects this |
| Feeding replicates to cosinor as one pooled single-fit | Use population-mean cosinor (subject = replication unit) so between-subject variance enters the CI and the test |

## The over-detection controversy (state it as live)

Laloum & Robinson-Rechavi (2020) showed that across seven popular methods (ARS, LS, RAIN, JTK, eJTK, GeneCycle, meta2d) rhythm calls are consistent and biologically meaningful ONLY for strong-amplitude signals; weak-signal calls are method-dependent and largely non-functional. There is no consensus "correct" method. The pragmatic (not full) response: (1) require an amplitude/rAMP effect-size filter IN ADDITION to FDR; (2) prefer methods with calibrated empirical nulls (eJTK, BooteJTK) over raw JTK on sparse data; (3) verify the genome-wide p-value histogram is roughly uniform with a spike near 0 before trusting any q. Report the amplitude distribution of the hit list, not just "N% of the transcriptome is rhythmic."

## Parameter Guide

| Parameter | Typical value | Rationale |
|-----------|---------------|-----------|
| Period | 24h (12h for ultradian) | Specified a priori; this is a test, not a search |
| Period window | 24 (LD) / 22-26 (DD) | Entrained locks to 24h; free-running tau != 24h. Wide windows inflate false positives |
| Sampling interval | 2-4h | Nyquist (<=12h) is a floor, not a target; shape resolution needs 1-2h |
| Cycles | >=2 (>=3 better) | One cycle cannot separate rhythm from trend/transient |
| Samples/cycle | >=6 (8-12+ better) | Six gives stable fit df; more resolves waveform and calibrates FDR |
| Replicates/timepoint | >=2-3 | Single replicate has no within-timepoint variance; FDR miscalibrates |
| FDR threshold | q < 0.05 | Necessary but not sufficient; always pair with an amplitude filter |
| Relative amplitude | rAMP > 0.1 | >=10% of baseline as a biological-relevance floor; a convention to sweep, not a law |

## Related Skills

temporal-genomics/differential-rhythmicity - Comparing rhythms between conditions (gain/loss/phase/amplitude change)
temporal-genomics/periodicity-detection - Unknown-period discovery with Lomb-Scargle and wavelets
temporal-genomics/temporal-clustering - Group rhythmic genes by phase/shape
differential-expression/timeseries-de - Temporal differential expression (a monotone trend, not rhythmicity)
data-visualization/heatmaps-clustering - Circular phase heatmaps and phase-ordered maps

## References

- Hughes ME, Hogenesch JB, Kornacker K. 2010. JTK_CYCLE: an efficient nonparametric algorithm for detecting rhythmic components in genome-scale data sets. J Biol Rhythms 25(5):372-380. doi:10.1177/0748730410379711
- Hughes ME, Abruzzi KC, Allada R, et al. 2017. Guidelines for genome-scale analysis of biological rhythms. J Biol Rhythms 32(5):380-393. doi:10.1177/0748730417728663
- Thaben PF, Westermark PO. 2014. Detecting rhythms in time series with RAIN. J Biol Rhythms 29(6):391-400. doi:10.1177/0748730414553029
- Wu G, Anafi RC, Hughes ME, Kornacker K, Hogenesch JB. 2016. MetaCycle: an integrated R package to evaluate periodicity in large scale data. Bioinformatics 32(21):3351-3353. doi:10.1093/bioinformatics/btw405
- Yang R, Su Z. 2010. Analyzing circadian expression data by harmonic regression based on autoregressive spectral estimation (ARSER). Bioinformatics 26(12):i168-i174. doi:10.1093/bioinformatics/btq189
- Hutchison AL, Maienschein-Cline M, Chiang AH, et al. 2015. Improved statistical methods enable greater sensitivity in rhythm detection for genome-wide data (eJTK). PLoS Comput Biol 11(3):e1004094. doi:10.1371/journal.pcbi.1004094
- Cornelissen G. 2014. Cosinor-based rhythmometry. Theor Biol Med Model 11:16. doi:10.1186/1742-4682-11-16
- Laloum D, Robinson-Rechavi M. 2020. Methods detecting rhythmic gene expression are biologically relevant only for strong signal. PLoS Comput Biol 16(3):e1007666. doi:10.1371/journal.pcbi.1007666
- Mei W, Jiang Z, Chen Y, Chen L, Sancar A, Jiang Y. 2021. Genome-wide circadian rhythm detection methods: systematic evaluations and practical guidelines. Brief Bioinform 22(3):bbaa135. doi:10.1093/bib/bbaa135
- Moškon M. 2020. CosinorPy: a python package for cosinor-based rhythmometry. BMC Bioinformatics 21:485. doi:10.1186/s12859-020-03830-w
