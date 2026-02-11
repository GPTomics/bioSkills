# Circadian Rhythm Detection - Usage Guide

## Overview

Detects circadian and ultradian rhythms in time-series omics data. Fits cosinor regression models to estimate rhythm parameters (amplitude, phase, MESOR) and applies non-parametric tests (JTK_CYCLE, RAIN) to identify genes with significant periodic expression at a specified period. Supports both Python (CosinorPy) and R (MetaCycle, RAIN, DiscoRhythm) workflows.

## Prerequisites

### Python
```bash
pip install cosinorpy pandas numpy statsmodels matplotlib
```

### R
```r
install.packages('data.table')
BiocManager::install(c('MetaCycle', 'rain', 'DiscoRhythm'))
```

### Data Requirements
- Expression matrix with genes as rows and timepoints as columns
- At least 2 complete cycles of the target period (e.g., 48h for circadian)
- Sampling interval no greater than half the target period (Nyquist criterion)
- Recommended: 6+ timepoints per cycle with biological replicates

## Quick Start

Tell the AI agent what to analyze:
- "Test which genes are circadian in my 48-hour time-course RNA-seq data"
- "Run JTK_CYCLE on my expression matrix to find rhythmic genes"
- "Fit cosinor models to my time-series data and extract phase and amplitude"
- "Compare circadian rhythmicity between control and treated conditions"

## Example Prompts

### Basic Rhythm Detection
> "I have a gene expression matrix sampled every 4 hours over 48 hours. Test which genes have circadian rhythms."

> "Run MetaCycle on my RNA-seq time-course data to identify 24-hour periodic genes."

### Parameter Estimation
> "Fit cosinor models to my circadian data and give me the phase, amplitude, and MESOR for each gene."

> "I need to know when each rhythmic gene peaks -- extract acrophase from my time-series expression data."

### Multi-Method Comparison
> "Compare JTK_CYCLE, RAIN, and cosinor results on my circadian dataset. Which genes are consistently rhythmic?"

> "Run both CosinorPy and MetaCycle on my data and find the overlap of significant rhythmic genes."

### Condition Comparison
> "I have time-course data for WT and KO mice. Which genes lose their circadian rhythm in the knockout?"

> "Compare the amplitude and phase of rhythmic genes between fed and fasted conditions."

## What the Agent Will Do

1. Load and validate the time-series expression matrix
2. Verify sufficient timepoints and cycles for the target period
3. Apply rhythmicity tests (cosinor, JTK_CYCLE, RAIN, or MetaCycle meta-analysis)
4. Correct p-values for multiple testing (BH FDR)
5. Extract rhythm parameters: amplitude, acrophase, period, MESOR
6. Filter significant rhythmic genes (q < 0.05)
7. Generate phase distribution plots and rhythmic gene heatmaps
8. Optionally compare rhythmicity between conditions

## Tips

- Two complete cycles (48h for circadian) is the absolute minimum; three cycles improves statistical power substantially
- JTK_CYCLE is the most widely used method for evenly sampled circadian data and is fast even genome-wide
- RAIN handles asymmetric waveforms (rapid induction, slow decay) better than cosinor
- Use MetaCycle meta2d() to integrate multiple methods and increase confidence in results
- Relative amplitude (rAMP) is more comparable across genes than raw amplitude since it normalizes by baseline expression
- Phase estimates are unreliable for genes near the significance threshold; focus interpretation on strongly rhythmic genes
- For unevenly sampled data, use Lomb-Scargle (see temporal-genomics/periodicity-detection) rather than JTK_CYCLE

## Related Skills

- temporal-genomics/periodicity-detection - Unknown-period discovery with Lomb-Scargle and wavelets
- temporal-genomics/temporal-clustering - Group rhythmic genes by phase
- differential-expression/timeseries-de - Temporal DE testing
- data-visualization/heatmaps-clustering - Circular phase heatmaps
