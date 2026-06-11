# MOFA2 Integration - Usage Guide

## Overview

MOFA2 is unsupervised: it generalizes PCA to multiple omics, decomposing variation into latent factors with per-view loadings and a variance-explained table (R^2 per factor per view). That decomposition is the deliverable. Because the factors are found without the phenotype, a factor need not align with the groups of interest - which makes "a factor that separates my cases and controls" meaningful evidence, but only if the factor is annotated, checked against technical covariates, and shown to be robust rather than mined from many candidates.

The skill owns the MOFA model itself - configuration, training, the variance decomposition, factor labeling, and signed-weight enrichment. Per-omic normalization and HVG selection happen upstream in data-harmonization; supervised integration (where the outcome drives the projection) is the sibling mixomics-analysis; the single-cell MOFA+ path lives in single-cell/multimodal-integration.

## Prerequisites

```r
BiocManager::install('MOFA2')
```
```bash
pip install mofapy2   # the Python training backend MOFA2 calls via basilisk/reticulate
```

Conceptual prerequisites and notes:
- Each view must be features-by-samples (the transpose of the usual layout), per-omic normalized and variance-stabilized, and HVG-filtered so feature counts are within an order of magnitude across views.
- Prefer transforming counts to a Gaussian likelihood over using the Poisson likelihood; raw counts in a Gaussian view make factor 1 the library-size factor.
- MOFA tolerates missing samples in a view, so a mosaic cohort can be passed directly rather than intersected to complete cases.
- The trained model serializes to an `.hdf5` file whose schema is MOFA2-version-coupled; R and Python (muon/mofax) read the same file but have different option defaults.

## Quick Start

Tell your AI agent what you want to do:
- "Integrate my RNA-seq, proteomics, and methylation to find shared factors"
- "Which factors explain the most variance, and in which views?"
- "Tell me which of my MOFA factors are batch factors, not biology"
- "Run enrichment on the weights of factor 2, separately for each sign"
- "Is factor 3 robust, or does it disappear when I change the seed?"

## Example Prompts

### Unsupervised factor discovery
> "I have RNA-seq, proteomics, and methylation on the same 60 patients, with methylation missing for 12 of them. Find the joint latent factors with MOFA2, tell me how much variance each factor explains per view, and which factors are shared versus view-specific."

### Labeling and excluding technical factors
> "Correlate my MOFA factors with both my clinical condition and my technical covariates (batch, sequencing depth), and tell me which factors I should exclude from biological interpretation because they track a technical variable."

### Signed-weight interpretation
> "Extract the top positive and top negative features for factor 1 in the RNA view, and run enrichment separately for each sign so I do not conflate the two poles of the axis."

### Robustness
> "Retrain my MOFA model with a different seed and fewer factors, and tell me whether the factor that separated my groups recurs - I do not want to build a story on a fragile factor."

## What the Agent Will Do

1. Prepare each view: confirm features-by-samples orientation, upstream normalization, and comparable feature counts.
2. Create the MOFA object and inspect the data-overview plot for the missing-data pattern.
3. Configure the model: over-specify the factor count and let ARD prune, match likelihoods (preferring Gaussian after transform), set a seed.
4. Train by variational inference and read the variance-explained table first to separate shared from view-specific factors.
5. Attach metadata after fitting and correlate every factor with biological and technical covariates, excluding batch/depth factors from interpretation.
6. Annotate the surviving factors with signed-weight enrichment and confirm they recur across seeds before reporting them.

## Key Concepts

- Views: the omics modalities (RNA, protein, methylation), each features-by-samples.
- Factors: unsupervised latent axes of variation; unordered and ARD-pruned, so factor 1 is not "most important" the way PC1 is.
- Weights: signed feature loadings; large positive and large negative weights are opposite poles of the same axis.
- Variance explained: R^2 per factor per view - the central output that says which factors are shared.
- Groups (MOFA+): a sample partition that lets factor activity differ across groups while weights stay shared. This is not batch correction and does not make MOFA supervised.

## Tips

- Read the variance-explained table before interpreting any factor; factors with near-zero R^2 everywhere are noise.
- Always correlate factors with technical covariates; a factor that tracks batch is a batch factor, not confounded biology.
- MOFA never saw the groups, so a factor that splits them is strong evidence - unless many factors were scanned to find it. Pre-specify or correct for the number of factors, and validate out-of-sample.
- Transform counts to Gaussian upstream rather than relying on the Poisson likelihood.
- Equalize feature counts across views by filtering the wider views; reach for `scale_views=TRUE` only when filtering cannot, and remember it changes the R^2 reading to within-view relative variance.
- Factor analysis needs sample size (the package floor is above 15); at tens of samples request fewer factors and validate out-of-sample.
- Confirm any headline factor recurs across seeds (factor-score correlation near 1); stochastic inference in particular is seed-dependent.

## Related Skills

- integration-design - The method-selection decision; MOFA is the default once correspondence is vertical
- mixomics-analysis - Supervised DIABLO/sPLS where the outcome drives the projection
- data-harmonization - Per-omic transform, HVG selection, and batch regression before MOFA
- similarity-network - Hard patient stratification alternative to soft factors
- single-cell/multimodal-integration - Single-cell MOFA+ (CITE-seq/Multiome)
- pathway-analysis/gsea - Enrichment of factor weights
- clinical-biostatistics/survival-analysis - Survival validation using factors as features
- workflows/multi-omics-pipeline - End-to-end multi-omics integration pipeline
