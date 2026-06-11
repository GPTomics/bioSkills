# Data Harmonization - Usage Guide

## Overview

Harmonization is the seam between the per-omic categories and the integration-method skills. It takes one already-normalized matrix per omic and emits a comparable, sample-matched set of blocks that MOFA2, mixOmics, or SNF can consume. It exists because a shared-latent integrator is blind to what an omic is: it sees one stacked matrix and spends its first factors on whichever block has the largest raw variance, so the transform, scaling, batch, and missing-value choices made here silently determine the integration result.

This skill owns the CROSS-omic decisions (which transform per omic, per-view versus per-feature scaling, the batch strategy, missing-data triage, the container). The DEEP normalization of any single omic (DESeq2 VST internals, methylation noob/BMIQ, proteomics VSN, metabolomics drift correction, compositional CLR theory) lives in that omic's own category and is cross-referenced, not re-derived here.

## Prerequisites

```r
BiocManager::install(c('MultiAssayExperiment', 'SummarizedExperiment', 'sva', 'limma'))
install.packages('imputeLCMD')   # MNAR-aware proteomics/metabolomics imputation
```

Conceptual prerequisites and notes:
- Each omic must already be per-omic normalized in its own category before harmonization (raw counts, raw beta values, or un-normalized intensities should not reach this skill).
- The blocks must share a sample axis; the MultiAssayExperiment sampleMap makes that linkage a structural invariant rather than string-matching.
- Two corrections happen in order: the per-omic variance-stabilizing transform (makes each block roughly Gaussian and homoscedastic) THEN per-view scaling (equalizes block contribution). Scaling a heteroscedastic block does not make it homoscedastic.
- Matrix orientation differs downstream: Bioconductor is samples-in-columns, while mixOmics and AnnData/MuData are samples-in-rows. Transpose on every cross-language hop.

## Quick Start

Tell your AI agent what you want to do:
- "Get my RNA, protein, and methylation blocks onto a common footing for MOFA2"
- "Which transform should each of my omics get before integration?"
- "Equalize my blocks so the 850k-CpG methylation does not dominate the factors"
- "Correct batch across my omics without deleting biology"
- "My proteomics is 30% missing below detection - how should I impute it?"

## Example Prompts

### Transform and scaling
> "I have RNA-seq counts, proteomics intensities, and 450k methylation on the same 60 samples. Apply the right transform to each block, then scale them so no single omic dominates the shared factors, and tell me which block was at risk of dominating."

### Batch strategy
> "My omics were generated in two batches. Should I run ComBat per omic before integrating, or model batch as a covariate inside the integration? First check whether batch is confounded with my case/control variable."

### Missing-data triage
> "My proteomics block has values missing below the detection limit and some sporadic gaps. Impute them by the right mechanism, and decide what to do about the 12 samples that have no proteomics at all."

### Container and matching
> "Assemble my three omics into a MultiAssayExperiment, tell me how many samples have all three assays, and decide whether to intersect to complete cases or keep the mosaic structure for MOFA2."

## What the Agent Will Do

1. Assemble the blocks into a MultiAssayExperiment and quantify the mosaic structure (how many samples have every omic).
2. Choose a variance-stabilizing transform per omic (or confirm one was applied upstream), routing the deep normalization to the per-omic category.
3. Equalize block contribution with per-view scaling, filtering near-constant features first so per-feature scaling cannot inflate noise.
4. Decide the batch strategy: cross-tabulate batch against the biological variable, refuse to correct a confounded design, and otherwise correct once per omic or model batch as a covariate - never both, never on a stacked matrix.
5. Triage missingness by mechanism (MNAR below detection vs MAR sporadic vs whole missing samples), avoiding whole-sample imputation in favor of a missing-tolerant integrator.
6. Emit a harmonized, sample-matched list of blocks for mofa-integration, mixomics-analysis, or similarity-network.

## Tips

- Apply the per-omic transform first, then scale across blocks; the two corrections are sequential and both required.
- Prefer per-view scaling over per-feature scaling for cross-block equalization, and filter near-constant features before any per-feature scaling.
- Cross-tabulate batch against the biological variable before correcting; if any cell is empty the design is confounded and the correction will delete biology.
- Model batch as a covariate for inferential steps; reserve scrubbing (ComBat, removeBatchEffect output) for visualization or tools that cannot accept a covariate. Correct once, in one place.
- Match imputation to the mechanism: QRILC/MinProb for below-detection MNAR gaps, kNN/missForest only for sporadic MAR gaps.
- Do not impute a whole missing omic to satisfy a complete-case method; MOFA2 models missing-view samples natively, so route mosaic cohorts there.
- Join omics on stable identifiers (Ensembl, UniProt, RefMet), never gene symbols, and decide the collapse rule for one-to-many mappings explicitly.

## Related Skills

- integration-design - The method-selection decision this harmonization feeds
- mofa-integration - Consumes harmonized blocks; models missing-view samples natively
- mixomics-analysis - Consumes harmonized blocks; needs complete cases and is per-feature scaled
- similarity-network - Consumes harmonized blocks for patient stratification
- differential-expression/batch-correction - Single-omic RNA-seq batch mechanics
- methylation-analysis/array-preprocessing - Methylation beta/M-value normalization
- proteomics/proteomics-qc - Proteomics normalization and QC
- metabolomics/normalization-qc - Metabolomics normalization and scaling
- metagenomics/abundance-estimation - Compositional/CLR theory for compositional omics
