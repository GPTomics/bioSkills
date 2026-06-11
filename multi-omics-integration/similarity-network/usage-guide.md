# Similarity Network Fusion - Usage Guide

## Overview

SNF integrates omics in patient-similarity space: it builds one patient-by-patient similarity network per omic and fuses them by cross-network diffusion, so feature count buys no votes and the n-much-smaller-than-p curse is dissolved into an n-by-n object. The fused graph is then spectral-clustered into candidate subtypes. The central caution is that spectral clustering always returns the requested number of clusters, so a subtype is a claim to defend with stability, survival separation, and replication - not something the algorithm discovers.

The skill owns SNF and the integrative-clustering landscape (NEMO, PINS, iCluster, CIMLR, intNMF, consensus) plus the validation discipline. Feature-space integration (which molecular axes or features matter) is the sibling mofa-integration or mixomics-analysis; survival mechanics live in clinical-biostatistics/survival-analysis; single-cell graph clustering is a different object in single-cell/clustering.

## Prerequisites

```r
install.packages(c('SNFtool', 'igraph', 'pheatmap'))
```

Conceptual prerequisites and notes:
- SNF requires complete data: every patient must have every omic, because the cross-diffusion multiplies aligned n-by-n matrices. For a mosaic cohort, use NEMO instead of dropping patients.
- Standardize each continuous omic per feature before computing distances, or a high-variance feature dominates the distance.
- The SNFtool API is easy to misread: `dist2` returns SQUARED Euclidean distance (take the square root), `affinityMatrix`'s width argument is `sigma` (not alpha), and `spectralClustering`'s `K` is the number of clusters, colliding with the K-neighbors argument elsewhere.
- SNF produces clusters but no feature-level model; feature attribution is post-hoc via `rankFeaturesByNMI`.

## Quick Start

Tell your AI agent what you want to do:
- "Identify patient subtypes by fusing my RNA-seq and methylation data"
- "How many subtypes does my fused network actually support?"
- "Show me whether the fusion beat my best single omic"
- "Check that my subtypes are stable when I resample patients"
- "Some patients are missing an omic - what should I use instead of SNF?"

## Example Prompts

### Patient stratification with defense
> "Fuse my RNA, protein, and methylation networks with SNF and propose patient subtypes. Then defend the subtype count: show me the eigengap estimates, whether the clusters are stable when I resample patients, and whether the fusion actually separates outcome better than my best single omic."

### Cluster-number honesty
> "Estimate the number of clusters from my fused network, but treat that as a plausibility check, not a truth. Show me how the clustering changes if I move K from 20 to 25, and tell me whether my headline subtypes survive."

### Partial data
> "About a third of my patients are missing proteomics. Should I still use SNF, or switch to a method built for partial data, and what does the intersection cost me if I force complete cases?"

### Fused versus single omic
> "Quantify whether my fused clustering agrees more with one omic than the others, and tell me honestly whether integrating helped or whether one omic was carrying the signal."

## What the Agent Will Do

1. Confirm complete multi-omic data on matched patients; for mosaic cohorts, recommend NEMO and report what the intersection would drop.
2. Standardize each omic, compute square-rooted distances, and build local-scaled affinity networks.
3. Fuse the networks with SNF and read the four cluster-number estimates as plausibility, not truth.
4. Spectral-cluster the fused graph into candidate subtypes at a defended cluster count.
5. Defend the subtypes: stability under resampling, a fused-versus-best-single-omic comparison, and covariate-adjusted survival separation.
6. Attribute features post-hoc with `rankFeaturesByNMI` and route survival mechanics and any biomarker claim out for validation.

## Key Parameters

| Parameter | SNFtool argument | Typical range |
|-----------|------------------|---------------|
| K (neighbors) | `affinityMatrix(diff, K, sigma)`, `SNF(Wall, K, t)` | 10-30 |
| sigma (kernel width) | `affinityMatrix(diff, K, sigma)` | 0.3-0.8 |
| t (fusion iterations) | `SNF(Wall, K, t)` | 10-20 |
| C (clusters) | `spectralClustering(W, C)` | the claim, defended |

## Tips

- Standardize each omic per feature before distances; `dist2` returns squared distance, so take the square root.
- Pass the kernel width as `sigma`, and remember `spectralClustering`'s K is the cluster count, not the neighbor count.
- Never tune K and sigma to maximize agreement with known labels and then claim that agreement; fix or pre-register a small grid and report stability across it.
- SNF needs complete data; switch to NEMO for mosaic cohorts rather than dropping or imputing whole omics.
- Always check whether the fusion beat the best single omic; multi-omics is not automatically better.
- Treat `rankFeaturesByNMI` output as post-hoc characterization, not a model that generated the clusters.
- Validate subtypes with covariate-adjusted survival and a permutation p, and replicate in an independent cohort before claiming biological subtypes.

## Related Skills

- integration-design - The method-selection and paired-vs-mosaic decision
- mofa-integration - Feature-space latent factors (which molecular axes matter)
- mixomics-analysis - Supervised feature signatures for known classes
- data-harmonization - Per-omic scaling before building distances
- clinical-biostatistics/survival-analysis - Survival validation of discovered subtypes
- single-cell/clustering - Graph clustering of cells (different object)
- pathway-analysis/go-enrichment - Enrichment of subtype-associated features
- workflows/multi-omics-pipeline - End-to-end multi-omics integration pipeline
