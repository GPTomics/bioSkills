# temporal-genomics

## Overview

Analyze temporal patterns in omics time-series data. Covers circadian rhythm detection, temporal gene clustering, trajectory modeling with GAMs, dynamic gene regulatory network inference, and periodicity discovery.

**Tool type:** mixed | **Primary tools:** CosinorPy, Mfuzz, mgcv, statsmodels, scipy

## Skills
| Skill | Description |
|-------|-------------|
| circadian-rhythms | Detect circadian and ultradian rhythms with CosinorPy and MetaCycle |
| temporal-clustering | Cluster genes by temporal expression profile with Mfuzz and DEGreport |
| trajectory-modeling | Model continuous temporal trajectories with GAMs and changepoint detection |
| temporal-grn | Infer dynamic gene regulatory networks from bulk time-series data |
| periodicity-detection | Discover periodic signals of unknown period with Lomb-Scargle and wavelets |

## Example Prompts
- "Test which genes have circadian expression patterns in my time-course data"
- "Cluster my temporally variable genes by expression profile shape"
- "Fit smooth curves to gene expression over time and compare conditions"
- "Infer time-delayed regulatory relationships between transcription factors and targets"
- "Find periodic patterns in my unevenly sampled time-series data"

## Requirements
```bash
# Python
pip install cosinorpy scipy pywt tslearn ruptures statsmodels

# R
install.packages(c('MetaCycle', 'mgcv', 'segmented', 'bnlearn'))
BiocManager::install(c('Mfuzz', 'DEGreport', 'rain', 'TCseq', 'DiscoRhythm'))

# dynGENIE3 (GitHub only)
devtools::install_github('vahuynh/dynGENIE3/dynGENIE3R')
```

## Related Skills

- **differential-expression** - Temporal DE testing with limma/maSigPro
- **gene-regulatory-networks** - Static and condition-specific GRN inference
- **single-cell** - Pseudotime trajectory analysis for single-cell data
- **pathway-analysis** - Functional enrichment of temporal gene clusters
- **data-visualization** - Plotting temporal profiles and heatmaps
