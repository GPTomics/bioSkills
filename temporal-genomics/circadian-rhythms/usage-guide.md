# Known-Period Rhythm Testing - Usage Guide

## Overview

Tests and estimates rhythmicity at a PRE-SPECIFIED period (canonically 24h) in a single condition of time-series omics. This is known-period TESTING (the period is fixed, then rhythmicity, phase, amplitude, and MESOR are tested/estimated at it), which is categorically different from unknown-period DISCOVERY (finding the period; see temporal-genomics/periodicity-detection) and from comparing rhythms BETWEEN conditions (differential rhythmicity; see temporal-genomics/differential-rhythmicity). It fits cosinor regression (CosinorPy), runs JTK_CYCLE/ARSER/Lomb-Scargle meta-analysis (MetaCycle meta2d), and non-parametric tests for asymmetric waveforms (RAIN, DiscoRhythm). The governing rule: temporal conclusions are dominated by SAMPLING DESIGN, not the algorithm.

## Prerequisites

### Python
```bash
pip install "numpy<2.0" CosinorPy pandas statsmodels matplotlib seaborn
```
CosinorPy 3.1 calls the removed `np.round_`, so it requires numpy<2.0, and it imports as `from CosinorPy import cosinor, cosinor1, file_parser` (capitalized package name).

### R
```r
install.packages(c('data.table', 'MetaCycle'))
BiocManager::install(c('rain', 'DiscoRhythm'))
```

### Data Requirements
- Feature x timepoint matrix (long format for CosinorPy: columns x=time, y=value, test=feature id)
- At least 2 complete cycles of the target period (48h for circadian; 3 cycles is better)
- >=6, ideally 8-12+, samples per cycle (Nyquist's 2/cycle is a floor, not a target)
- >=2-3 biological replicates per timepoint for calibrated FDR
- A declared light regime (LD entrained -> ZT; DD free-running -> CT) and period window

## Quick Start

Tell the AI agent what to analyze:
- "Test which genes are circadian in my 48-hour time-course RNA-seq (sampled every 4 hours)"
- "Fit cosinor models and give me phase, amplitude, MESOR, and relative amplitude per gene"
- "Run MetaCycle meta2d on my expression matrix and rank rhythmic genes with an amplitude filter"
- "My waveforms are asymmetric (fast rise, slow decay) - use RAIN instead of cosinor"

## Example Prompts

### Basic Rhythm Detection
> "I have a gene expression matrix sampled every 4 hours over 48 hours. Test which genes have 24-hour rhythms and report relative amplitude, not just p-value."

> "Run MetaCycle meta2d on my RNA-seq time-course, but treat meta2d_BH.Q as a rank and add an rAMP filter."

### Parameter Estimation
> "Fit cosinor models and give me acrophase converted to peak hour, amplitude, and MESOR for each rhythmic gene."

> "This is free-running (constant darkness) data - use a 22-26h period window and report phase in CT."

### Method Choice
> "My clock-output genes have sharp asymmetric peaks. Which method should I use and why?"

> "I have single-replicate sparse sampling - are my JTK p-values trustworthy, and what should I use instead?"

(Comparing rhythms between conditions - "which genes lose/gain/phase-shift their rhythm in the KO?" - is temporal-genomics/differential-rhythmicity, not this skill.)

## What the Agent Will Do

1. Confirm the question is known-period testing (not period discovery) and set the period window from the light regime
2. Load and validate the matrix; check cycles, sampling density, and replication against design minima
3. Choose a method by waveform shape and sampling (cosinor / JTK / eJTK / ARSER / RAIN / meta2d) and explain why
4. Fit models and extract amplitude, relative amplitude, acrophase (converted to peak hour), MESOR, p-value
5. Control FDR (BH), then apply an effect-size filter (rAMP / fold-change) because significance alone over-detects
6. Report the amplitude distribution of the hit list and flag design confounds (harvest order, LD masking, damping)
7. Route between-condition questions to temporal-genomics/differential-rhythmicity (never intersect two separate rhythm lists)

## Tips

- Fix the period window to the light regime: 24h for entrained (LD) data, ~22-26h for free-running (DD, where tau != 24h). A wide window turns a test into a discovery search and inflates false positives.
- CosinorPy's fit_group already returns a BH-adjusted q column; recompute BH only to control the correction set, and pass method='fdr_bh' explicitly (the statsmodels default is Holm-Sidak).
- Convert CosinorPy acrophase to peak hour with (-acrophase)*T/(2*pi) mod T, and sanity-check against a known clock gene (mouse liver Arntl peaks ~CT22-0, Dbp ~CT8-10).
- Relative amplitude (rAMP = amplitude/MESOR) is comparable across genes; raw amplitude scales with expression level and is not.
- Significance is necessary but not sufficient (Laloum 2020): add an rAMP or fold-change filter and report the amplitude distribution, not just "N% of the transcriptome is rhythmic."
- On sparse single-replicate designs, JTK p-values are quantized and anti-conservative; prefer eJTK/BooteJTK and inspect the genome-wide p-value histogram before trusting FDR.
- meta2d_BH.Q is Fisher over correlated ARS/JTK/LS p-values, so read it as a consensus rank, not a literal FDR; distrust the averaged period/phase when the constituent methods disagree.
- RAIN detects asymmetric waveforms (fast induction, slow decay) that cosinor and JTK miss, but has lower power for symmetric sinusoids.
- For between-condition comparison, use temporal-genomics/differential-rhythmicity (model-based LimoRhyde/dryR/compareRhythms) - never define "lost rhythm" as the set difference of two thresholded single-condition lists.
- Reduced BULK amplitude can mean cell desynchrony, not arrhythmia; report "reduced ensemble amplitude" and use single-cell or imaging assays to separate the causes.
- An apparent 24h rhythm under LD may be light/feeding-driven (masking); only persistence in DD demonstrates an endogenous circadian rhythm. Diurnal != circadian.
- Randomize sample processing order: harvest-order drift aliases perfectly onto circadian time and cannot be removed analytically.

## Related Skills

temporal-genomics/differential-rhythmicity - Comparing rhythms between conditions (gain/loss/phase/amplitude change)
temporal-genomics/periodicity-detection - Unknown-period discovery with Lomb-Scargle and wavelets
temporal-genomics/temporal-clustering - Group rhythmic genes by phase/shape
differential-expression/timeseries-de - Temporal differential expression (monotone trend, not rhythmicity)
data-visualization/heatmaps-clustering - Circular phase heatmaps and phase-ordered maps
