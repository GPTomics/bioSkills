# Multi-Omics Integration Design - Usage Guide

## Overview

This skill is the decision layer that runs BEFORE any integration tool. It maps a biological question to a method class, names the sample correspondence (vertical/paired, mosaic, horizontal, or diagonal), enforces the small-n / large-p discipline that makes a held-out cohort the real endpoint, and runs the per-view variance-imbalance diagnostic that exposes a "shared" factor that is really one omic in disguise. It exists because the most consequential multi-omics mistakes are made in the choice of strategy, not in the syntax of a tool: feeding horizontal (same-feature, multi-cohort) data to a vertical method, reporting an unvalidated in-cohort signature as a biomarker, or mistaking the dominance of the biggest omic for integration.

It is method-agnostic and routes to the tool skills: mofa-integration (unsupervised shared factors), mixomics-analysis (supervised DIABLO signatures, sPLS pairs, MINT multi-study), and similarity-network (patient stratification by network fusion). Single-cell multimodal integration is a different paradigm owned by single-cell/multimodal-integration.

## Prerequisites

```r
BiocManager::install(c('MultiAssayExperiment', 'SummarizedExperiment'))
```

Conceptual prerequisites and notes:
- The omics blocks must be measured on a SHARED sample axis (the same patients/specimens). If instead the same features were measured in different cohorts, that is horizontal integration (meta-analysis / batch correction), not this category.
- Per-omic normalization happens FIRST, in each omic's own category (RNA-seq -> VST/logCPM; methylation -> M-values; proteomics -> log2 + MNAR-aware imputation; metabolomics -> log + Pareto). This skill assumes each block is already on a sane scale; the cross-omic harmonization decision is owned by data-harmonization.
- Bulk multi-omics is almost always small-n / large-p: tens of samples and 10^4-10^6 features. That regime, not the choice of tool, dominates what is trustworthy.
- A MultiAssayExperiment is the natural container for the paired-vs-mosaic decision; matrix orientation differs by downstream tool (MOFA2 wants features x samples, mixOmics wants samples x features).

## Quick Start

Tell your AI agent what you want to do:
- "I have RNA, protein, and methylation on the same patients - which integration method fits?"
- "Is this vertical or horizontal integration?"
- "My cohort is mosaic - should I intersect to complete cases or model the missingness?"
- "Check whether one omic dominates my integration before I trust the factors"
- "How do I validate an integrated subtype - is in-cohort cross-validation enough?"

## Example Prompts

### Method selection from the question
> "I have transcriptomics, proteomics, and methylation on 70 tumors and want to discover molecular subtypes with no predefined labels. Which method class should I use, and what is the honest way to claim the subtypes are real rather than a clustering artifact?"

### Correspondence and mosaic decision
> "Three RNA-seq cohorts of the same disease, same genes - my collaborator wants to 'integrate' them with MOFA. Is that the right tool, or is this a different kind of integration? If different, where should it go?"

### The n<<p validation discipline
> "I have a DIABLO signature that classifies responders with 0.95 AUC in cross-validation on 44 samples. How much should I trust it, and what would make it credible?"

### Variance-imbalance diagnostic
> "Before I interpret my MOFA factors, check whether my 850k-CpG methylation block is dominating the shared latent space relative to my 200-metabolite block, and tell me what to do if it is."

## What the Agent Will Do

1. Establish the correspondence: confirm the blocks share a sample axis (vertical) versus share features across cohorts (horizontal), and assemble a MultiAssayExperiment to quantify how mosaic the cohort is.
2. Map the question to a method class: subtype discovery (SNF/iCluster), shared axes of variation (MOFA2/JIVE), predictive signature (DIABLO), or pairwise correlation (sPLS), routing single-cell and horizontal cases out.
3. Surface the n<<p regime: count samples versus features, and set the expectation that tuning needs repeated cross-validation and that a held-out cohort is the endpoint.
4. Run the variance-imbalance diagnostic: compare per-block variance and feature counts, and plan to read the post-fit per-view variance-explained table; recommend equalization or SNF if one omic dominates.
5. Hand off to the chosen tool skill with the correspondence, scaling, and validation constraints made explicit.

## Tips

- Name the correspondence before naming a tool. If the answer is "same features, different cohorts," stop - that is meta-analysis, not vertical integration.
- Treat an unvalidated integrated signature as hypothesis-generating, never as a biomarker. At small n the in-cohort fit is guaranteed and uninformative.
- Always read the per-view variance-explained table. If every shared factor loads on one omic, the integration re-discovered the biggest data type; equalize the blocks or move to SNF.
- Prefer modeling the missingness (MOFA2) over intersecting a mosaic cohort to complete cases, because discarding incomplete samples is expensive when n is already small.
- There is no universally best integration method; benchmarks disagree by design. Cross-check a headline result with a second method class - a real subtype should appear in both a MOFA map and an SNF clustering.
- Choose the DIABLO design matrix from the question a priori (discrimination-first low, mechanistic-correlation-first high) and report the classification cost; tuning it until the result looks "integrated" is p-hacking by hyperparameter.

## Related Skills

- mofa-integration - Unsupervised shared-factor discovery (the default tool once correspondence is vertical)
- mixomics-analysis - Supervised DIABLO signatures, sPLS pairs, and MINT multi-study integration
- similarity-network - Patient stratification via similarity network fusion
- data-harmonization - Per-block normalization, scaling, batch, and the mosaic missing-omic decision
- single-cell/multimodal-integration - Single-cell CITE-seq/Multiome integration (different paradigm)
- differential-expression/batch-correction - Horizontal same-feature meta-analysis and batch correction
- machine-learning/model-validation - Cross-validation and overfitting theory for supervised integration
- clinical-biostatistics/survival-analysis - Survival validation of discovered subtypes
- workflows/multi-omics-pipeline - End-to-end multi-omics integration pipeline
