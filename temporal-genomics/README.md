# temporal-genomics

## Overview

Analyze temporal patterns in omics time-series data, with the postdoc-level guardrails that make the calls trustworthy: sampling-design constraints (Nyquist floor, >=2 cycles, evenness), the multiple-testing and null-model pitfalls specific to rhythm/periodicity screens, and the causal-inference ceiling of bulk temporal networks. Covers known-period rhythm detection and differential rhythmicity, unknown-period discovery, temporal profile clustering, GAM/changepoint trajectory modeling, and dynamic gene regulatory network inference.

**Tool type:** mixed | **Primary tools:** CosinorPy, MetaCycle, scipy/astropy, Mfuzz, mgcv, statsmodels

## Skills
| Skill | Description |
|-------|-------------|
| circadian-rhythms | Test known-period (usually 24h) rhythmicity and estimate phase/amplitude in a single condition with cosinor, JTK_CYCLE/eJTK, ARSER, RAIN, MetaCycle; rAMP effect-size filter against over-detection |
| differential-rhythmicity | Compare rhythms BETWEEN conditions/genotypes/tissues (gain/loss/phase/amplitude change) with LimoRhyde, dryR, compareRhythms; the detect-then-Venn anti-pattern and the DE-vs-DR distinction |
| periodicity-detection | Discover unknown-period signals with Lomb-Scargle (generalized, uneven-sampling), wavelets (cone-of-influence, red-noise significance), autocorrelation, and Welch PSD |
| temporal-clustering | Group genes by temporal profile shape with Mfuzz soft c-means, TCseq, DEGreport, and DTW (tslearn) after selecting temporally variable genes |
| trajectory-modeling | Model continuous bulk trajectories with mgcv GAMs (NB for counts, autocorrelation-aware) and detect abrupt regime shifts with segmented/ruptures changepoints |
| temporal-grn | Infer time-delayed, hypothesis-grade regulatory edges from bulk time series with Granger causality, dynGENIE3, and dynamic Bayesian networks |

## Example Prompts
- "Test which genes have circadian expression patterns in my time-course data"
- "Test whether rhythms differ between my wild-type and knockout mice"
- "Find periodic patterns of unknown period in my unevenly sampled time-series data"
- "Cluster my temporally variable genes by expression profile shape"
- "Fit smooth curves to gene expression over time and compare conditions"
- "Infer time-delayed regulatory relationships between transcription factors and targets"

## Requirements
```bash
# Python
pip install cosinorpy scipy astropy pywavelets tslearn ruptures statsmodels scikit-learn
# Note: CosinorPy 3.1 requires numpy<2.0 (it calls the removed np.round_)

# R (CRAN)
install.packages(c('MetaCycle', 'mgcv', 'segmented', 'bnlearn', 'limorhyde'))

# R (Bioconductor)
BiocManager::install(c('Mfuzz', 'DEGreport', 'rain', 'TCseq', 'DiscoRhythm', 'tradeSeq'))

# GitHub only
devtools::install_github('vahuynh/dynGENIE3/dynGENIE3R')     # dynamic GRN
devtools::install_github('naef-lab/dryR')                    # differential rhythmicity
devtools::install_github('bharathananth/compareRhythms')     # differential rhythmicity
```

## Related Skills

- **differential-expression** - Temporal DE testing with limma splines / DESeq2 LRT (select temporally variable genes before clustering)
- **gene-regulatory-networks** - Static and condition-specific GRN inference; single-cell regulons (SCENIC)
- **single-cell** - Pseudotime trajectory inference for single-cell data (distinct from bulk real-time trajectory modeling)
- **pathway-analysis** - Functional enrichment of temporal gene clusters (use the temporal gene set as background)
- **data-visualization** - Plotting temporal profiles, phase heatmaps, and scalograms
- **workflows/timecourse-pipeline** - End-to-end bulk time-course pipeline orchestrating these skills
