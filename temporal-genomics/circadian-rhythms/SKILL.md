---
name: bio-temporal-genomics-circadian-rhythms
description: Detects circadian and ultradian rhythms in time-series omics data using CosinorPy cosinor models, MetaCycle (JTK_CYCLE, ARSER), and RAIN non-parametric tests. Fits cosine models to estimate phase and amplitude, tests rhythmicity significance at pre-specified periods. Use when testing for 24-hour or other known-period oscillations in circadian, feeding-fasting, or light-dark cycle experiments. Not for unknown-period discovery (see temporal-genomics/periodicity-detection).
tool_type: mixed
primary_tool: CosinorPy
---

# Circadian Rhythm Detection

Identifies periodic gene expression patterns at known periods (typically 24h) using cosinor regression, non-parametric rhythmicity tests, and meta-analysis approaches combining multiple methods.

## Core Workflow

1. Prepare time-series expression matrix (genes x timepoints)
2. Fit cosinor models or apply rhythmicity tests at specified period
3. Extract rhythm parameters: amplitude, phase (acrophase), p-value
4. Correct for multiple testing (BH FDR)
5. Filter significant rhythmic genes and characterize phase distribution

## CosinorPy (Python)

### Single-Component Cosinor

Fits `y = M + A*cos(2*pi*t/T + phi)` where M = MESOR, A = amplitude, phi = acrophase.

```python
import pandas as pd
from cosinorpy import file_parser, cosinor, cosinor1

df = file_parser.read_csv('expression_timecourse.csv')

# Single-component cosinor fit for one gene
# period=24: standard circadian period in hours
# fit_me takes X (time) and Y (expression) arrays, returns a tuple
gene_data = df[df['test'] == 'test_gene']
res = cosinor.fit_me(gene_data['x'].values, gene_data['y'].values, period=24, n_components=1)
```

### Multi-Component Cosinor

Fits harmonics to capture non-sinusoidal waveforms (e.g., sharp peaks).

```python
# n_components=2: adds 12h harmonic to capture asymmetric waveforms
res_multi = cosinor.fit_me(gene_data['x'].values, gene_data['y'].values, period=24, n_components=2)

# n_components=3: adds 8h harmonic; rarely needed unless waveform is highly complex
res_3 = cosinor.fit_me(gene_data['x'].values, gene_data['y'].values, period=24, n_components=3)
```

### Population-Mean Cosinor

Combines individual fits across biological replicates or subjects.

```python
# Population-mean cosinor across multiple subjects
# Estimates group-level amplitude/phase with confidence intervals
pop_results = cosinor1.population_fit_cosinor(
    df, period=24, save_to='results/'
)
```

### Batch Processing

```python
# Genome-wide cosinor analysis via fit_group
# fit_group expects long-format DataFrame with columns: 'x' (time), 'y' (expression), 'test' (gene name)
results_df = cosinor.fit_group(df, period=24)

# BH correction for multiple testing (fit_group output has 'p' column; no built-in q-values)
from statsmodels.stats.multitest import multipletests
reject, qvals, _, _ = multipletests(results_df['p'], method='fdr_bh')
results_df['q_value'] = qvals
# q < 0.05: standard FDR threshold for rhythmicity
rhythmic = results_df[results_df['q_value'] < 0.05]
```

## MetaCycle (R)

Integrates JTK_CYCLE, ARSER, and Lomb-Scargle into a single meta-analysis framework.

```r
library(MetaCycle)

# meta2d integrates results from multiple methods
# cycMethod='JTK': JTK_CYCLE is the most widely used; robust and non-parametric
# minper=20, maxper=28: search window around 24h; allows detection of near-circadian periods
# timepoints: actual sampling times in hours (must match column order)
meta2d(
    infile = 'expression_matrix.csv',
    filestyle = 'csv',
    outdir = 'metacycle_results/',
    timepoints = seq(0, 44, by = 4),
    cycMethod = c('JTK', 'ARS', 'LS'),
    minper = 20,
    maxper = 28,
    outputFile = TRUE,
    outRawData = FALSE
)

# JTK_CYCLE alone (faster, suitable for evenly sampled data)
meta2d(
    infile = 'expression_matrix.csv',
    filestyle = 'csv',
    outdir = 'jtk_results/',
    timepoints = seq(0, 44, by = 4),
    cycMethod = 'JTK',
    minper = 24,
    maxper = 24
)
```

### MetaCycle Output Interpretation

```r
results <- read.csv('metacycle_results/meta2d_expression_matrix.csv')

# meta2d_pvalue: Fisher integration of per-method p-values
# meta2d_BH.Q: BH-corrected q-value across all genes
# meta2d_period: estimated period (hours)
# meta2d_phase: estimated peak time (hours from ZT0)
# meta2d_AMP: amplitude (half peak-to-trough distance)
# meta2d_rAMP: relative amplitude (AMP / baseline); robust across expression levels

# q < 0.05: standard FDR threshold
rhythmic <- results[results$meta2d_BH.Q < 0.05, ]
```

## RAIN (R/Bioconductor)

Non-parametric test that handles asymmetric waveforms (e.g., rapid induction, slow decay).

```r
library(rain)

# expression_mat: genes x timepoints matrix
# period=24: test for 24h rhythmicity
# deltat=4: sampling interval in hours
# nr.series=2: number of biological replicates per timepoint
# nr.series=2: samples per timepoint are grouped as replicates
results <- rain(
    t(expression_mat),
    period = 24,
    deltat = 4,
    nr.series = 2,
    method = 'independent'
)

# Adjust for multiple testing
results$q_value <- p.adjust(results$pVal, method = 'BH')
# q < 0.05: standard FDR threshold for rhythmicity detection
rhythmic <- results[results$q_value < 0.05, ]
```

## DiscoRhythm (R/Bioconductor)

Comprehensive framework with both scripted and interactive (Shiny) interfaces.

```r
library(DiscoRhythm)

# Scripted analysis pipeline
se <- discoGetSimu(TRUE)
disco_results <- discoBatch(se, report = NULL, osc_period = 24)
```

## Parameter Guide

| Parameter | Typical Value | Rationale |
|-----------|---------------|-----------|
| Period | 24h | Standard circadian period; use 12h for ultradian |
| Sampling interval | 2-4h | Nyquist: must sample at least 2x per period (every 12h minimum) |
| Minimum cycles | >=2 complete cycles | One cycle cannot distinguish trend from oscillation |
| Minimum timepoints | >=6 per cycle | JTK_CYCLE requires >=6; more improves power |
| FDR threshold | q < 0.05 | Standard multiple testing correction |
| Amplitude cutoff | Context-dependent | Often top 25th percentile of amplitudes among significant genes |
| Relative amplitude | rAMP > 0.1 | Filters low-amplitude oscillations; 10% of baseline is minimal biological relevance |

## Method Selection

| Method | Best For | Limitations |
|--------|----------|-------------|
| Cosinor | Parametric estimation, sinusoidal waveforms | Assumes sinusoidal shape |
| JTK_CYCLE | Robust non-parametric, evenly sampled | Requires even sampling |
| ARSER | Handles non-sinusoidal, spectral decomposition | Slower, needs longer series |
| RAIN | Asymmetric waveforms | Less power for symmetric waves |
| Lomb-Scargle | Uneven sampling | See periodicity-detection for unknown periods |

## Related Skills

- periodicity-detection - Unknown-period discovery with Lomb-Scargle and wavelets
- temporal-clustering - Group rhythmic genes by phase
- differential-expression/timeseries-de - Temporal DE testing
- data-visualization/heatmaps-clustering - Circular phase heatmaps
