# Batch Correction - Usage Guide

## Overview

Decision-grade batch correction for pooled CRISPR screens. Covers diagnosis (PCA + variance decomposition to decide whether batch dominates), the four primary correction strategies (ComBat empirical-Bayes, RUV, SVA, NTC-anchored normalization), the model-based alternative of including batch as a MAGeCK MLE / Chronos covariate (preferred for hit calling), screen-specific batch sources (cell-passage cohort, library lot, infection day, sequencing run, FBS lot, Cas9 lot), and the failure modes when correction destroys biology.

## Prerequisites

```bash
pip install pycombat pandas numpy scipy scikit-learn matplotlib seaborn
pip install mageck
# R (for SVA, RUVSeq, limma)
R -e "BiocManager::install(c('sva', 'RUVSeq', 'limma'))"
```

Required inputs: count matrix (rows=sgRNA, columns=samples), metadata table with at least `batch`, `condition`, and `replicate` columns, and (for NTC-anchored norm) a list of non-targeting control sgRNAs.

## Quick Start

Tell the AI agent what to do:
- "Diagnose batch effects in my screen: PCA, variance decomposition, is correction needed?"
- "Apply ComBat correction with biological_covariate=condition; verify PR-AUC is preserved post-correction"
- "Add batch as a covariate to MAGeCK MLE design matrix instead of pre-correcting"
- "Compare ComBat vs NTC-anchored normalization on my multi-batch screen"
- "Decide if I should correct or redesign because batch is fully confounded with condition"

## Example Prompts

### Diagnostics

> "Run PCA on my multi-batch screen colored by batch and condition. Compute F-statistics for batch and condition in each PC. Decide whether batch dominates (correction needed) or condition dominates (no correction needed)."

> "Replicate Pearson within batches is 0.95+ but across-batch is 0.78. Confirm batch effects are present; recommend ComBat vs NTC-anchored vs covariate modeling."

> "Check whether batch is confounded with condition. If yes (e.g., all drug arm in batch 2), correction will destroy biology -- explicit covariate modeling required instead."

### ComBat Application

> "Apply ComBat to log10(counts+1) with biological_covariate=condition. Verify post-correction PCA shows batches overlap and CEGv2 PR-AUC is preserved."

> "Compare ComBat with mod (biological covariate) vs without. Quantify how much biological signal is preserved in each case."

### Batch-Aware MLE

> "Build a MAGeCK MLE design matrix with explicit batch indicator columns. Run mageck mle on the multi-batch screen and report per-condition beta scores after batch adjustment. Compare to ComBat-then-MAGeCK approach."

> "Add SVA-discovered latent factors as covariates in MAGeCK MLE design matrix to capture unknown batch confounders."

### NTC-Anchored Normalization

> "Use the 800 non-targeting controls in my library as normalization anchor. Scale each sample so NTC median = 1000. Then run hit calling on NTC-anchored counts."

### Multi-Cell-Line Cancer Panels

> "I have 12 cancer cell lines screened across 4 batches. Switch to Chronos (built-in batch + CN modeling) instead of ComBat + MAGeCK; explain why this is preferred."

### Failure Diagnostics

> "After ComBat, my essentialome PR-AUC dropped from 0.78 to 0.41. Diagnose: did ComBat eliminate biological signal? Re-run with mod covariate or revert and use MLE-with-batch-covariate."

## What the Agent Will Do

1. Load count matrix and metadata; verify batch and condition columns
2. Run PCA + variance decomposition to diagnose batch dominance
3. Check confounding: is batch correlated with condition? If yes, recommend covariate modeling over pre-correction
4. Decide correction method based on diagnostic and batch annotation:
   - Annotated batches + ≥3 samples per batch -> ComBat with biological covariate
   - Unknown batch sources -> RUV with NTC reference set
   - Multi-screen with same library -> JACKS or Chronos (built-in batch handling)
   - Confounded batch+condition -> MLE with batch covariate (no pre-correction)
5. Apply correction or build MLE design matrix
6. Re-run PCA to confirm batches now overlap
7. Compute pre- vs post-correction CEGv2 PR-AUC (should be same or higher post)
8. Run hit calling on corrected counts; compare to uncorrected hit list
9. Output corrected counts, PCA before/after plots, correction summary

## Tips

- The single most common mistake is applying ComBat without biological covariate, which removes condition variance along with batch. Always supply `mod`.
- ComBat shifts counts before testing; the MLE-with-covariate approach correctly propagates batch-term uncertainty into the condition-beta standard error. ComBat-then-test is over-confident in FDR.
- Verify post-correction by re-checking the same metrics that motivated correction (PCA, PR-AUC, replicate Pearson). Don't trust correction blindly.
- For DepMap-style cancer cell-line panels, Chronos handles batch + CN bias + screen quality jointly. Use it instead of ComBat + MAGeCK.
- NTC-anchored normalization requires ≥500 NTCs in the library for stable median. Below this, fall back to median normalization.
- For in-vivo screens, batch sources are different (animal cohort, surgical day, tissue dissociation prep); see [[in-vivo-screens]] for in-vivo-specific batch handling.
- Sequential corrections (ComBat then RUV, or RUV then median) double-correct and destroy biology. Pick one method.
- Cancer-line screens have CN bias as a "batch-like" effect per cell line; CN correction (CRISPRcleanR / Chronos) is separate from batch correction. Apply CN correction first, then batch (or use Chronos which handles both).

## Decision Cheat Sheet

| Diagnostic finding | Action |
|--------------------|--------|
| PC1 batch F >> condition F | Apply ComBat with `mod` |
| Batch ⊥ condition | ComBat or MLE-with-covariate |
| Batch fully ⊆ condition | Redesign or only MLE-with-covariate |
| Multi-cell-line cancer | Chronos (preferred) |
| Multi-screen with same library | JACKS (joint efficacy) |
| Day-0 separates by batch but endpoint doesn't | Pre-screen technical bias; ComBat |
| Endpoint separates by batch but Day-0 doesn't | Selection-driven (FBS, drug lot); covariate model |
| <3 samples per batch | Cannot estimate batch; use covariate or redesign |

## Validation Checklist

After applying correction:
- [ ] PCA: batches now overlap (visual)
- [ ] Within-batch Pearson preserved (should be unchanged)
- [ ] Across-batch Pearson improved
- [ ] CEGv2 PR-AUC preserved or higher
- [ ] NTC distribution stable across batches
- [ ] No new outlier samples introduced

## Related Skills

- crispr-screens/mageck-analysis - Batch-aware MLE design matrix (preferred)
- crispr-screens/screen-qc - Pre-correction PCA + variance decomposition
- crispr-screens/copy-number-correction - Chronos handles batch + CN jointly
- crispr-screens/library-design - NTC composition for NTC-anchored norm
- crispr-screens/hit-calling - Post-correction hit calling
- crispr-screens/jacks-analysis - Joint analysis across batches
- crispr-screens/in-vivo-screens - In-vivo-specific batch sources
