# Periodicity Detection - Usage Guide

## Overview

Discovers periodic signals of unknown period in time-series omics data. Uses Lomb-Scargle periodograms for irregularly sampled data, autocorrelation for fundamental period estimation, and wavelet time-frequency decomposition for transient or time-varying periodicity. Designed for exploratory frequency analysis when the oscillation period is not known a priori.

## Prerequisites

### Python
```bash
pip install scipy numpy pywt astropy statsmodels matplotlib
```

### Data Requirements
- Time-series expression values with corresponding observation times
- At least 2 complete cycles of the shortest period of interest
- For Lomb-Scargle: observation times can be irregularly spaced
- For wavelet analysis: evenly spaced data preferred (interpolate if needed)

## Quick Start

Tell the AI agent what to search for:
- "Find periodic expression patterns in my unevenly sampled time-series data"
- "What are the dominant oscillation periods in my gene expression time course?"
- "Detect transient periodicity in my developmental time-series data"
- "Screen my RNA-seq time course for genes with significant periodic expression"

## Example Prompts

### Unknown Period Discovery
> "I have 72 hours of gene expression data sampled at irregular intervals. Find which genes have periodic expression and what their periods are."

> "Compute Lomb-Scargle periodograms for all genes in my time-course data and identify the dominant period for each."

### Cell Cycle Oscillations
> "Search for cell-cycle-related periodicities in my expression data. I expect periods around 18-24 hours but am not sure of the exact period."

> "Which genes oscillate with a period different from 24 hours in my time-series experiment?"

### Transient Periodicity
> "Some genes in my data seem to oscillate early but stop later. Use wavelet analysis to detect transient periodicity."

> "Apply continuous wavelet transform to my time-series expression data and show me a scalogram of the time-frequency structure."

### Genome-Wide Screening
> "Screen all 15,000 expressed genes for significant periodicity using Lomb-Scargle with FDR correction."

> "Run a permutation-based periodicity test on my time-course data and report genes with q < 0.05."

## What the Agent Will Do

1. Load time-series expression data and verify sampling structure
2. Detrend data if long-term trends are present (optional)
3. Compute frequency spectra (Lomb-Scargle, Welch, or autocorrelation)
4. Identify significant spectral peaks and estimate dominant periods
5. Assess statistical significance (false alarm probability or permutation FDR)
6. For transient periodicity, apply wavelet CWT and generate scalograms
7. Correct for multiple testing in genome-wide analyses
8. Export significant periodic genes with estimated periods and p-values

## Tips

- Lomb-Scargle is the default choice for biological time series because it handles irregular sampling natively
- Always check that the data span at least 2 complete cycles of the period of interest; one cycle cannot distinguish a trend from an oscillation
- Detrend (subtract polynomial fit) before frequency analysis if there is a strong monotonic trend; trends create low-frequency power that can mask true oscillations
- False alarm probability (FAP) from astropy is more reliable than naive significance thresholds on periodogram power
- Wavelet analysis trades frequency precision for time localization; use it when periodicity is expected to be transient
- For cell-cycle analysis, the expected period depends on cell type and growth conditions; do not assume 24h without evidence
- Permutation FDR is computationally expensive for genome-wide screening (1000 permutations x 15,000 genes); consider pre-filtering to variable genes first

## Related Skills

- temporal-genomics/circadian-rhythms - Known-period rhythm testing with cosinor and JTK_CYCLE
- temporal-genomics/temporal-clustering - Group genes by periodicity characteristics
- temporal-genomics/trajectory-modeling - Non-periodic trajectory fitting with GAMs
