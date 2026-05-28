# Differential Analysis - Usage Guide

## Overview
Differential analysis identifies cell populations that differ in abundance (DA) or marker expression (DS) between conditions. The defining principle this skill encodes is that the SAMPLE/subject - not the cell - is the experimental unit: diffcyt aggregates cells to per-sample-per-cluster counts (DA) and arcsinh-medians (DS), then tests across biological replicates with edgeR/limma/GLMM. It also covers the compositionality of cluster proportions (a real expansion forces apparent depletion elsewhere), mixed models for paired designs, modeling batch rather than cleaning it, and when to reach for simplex-aware compositional methods.

## Prerequisites
```bash
# R/Bioconductor
BiocManager::install(c('diffcyt', 'CATALYST', 'edgeR', 'limma'))
# Compositional alternatives (optional)
BiocManager::install(c('sccomp'))
```

## Quick Start
Tell your AI agent what you want to do:
- "Test which clusters differ in abundance between treatment and control"
- "Test Ki67 expression within each cluster between conditions"
- "Run a paired differential analysis for my pre/post samples"
- "Check whether my depletion result is a compositional artifact"

## Example Prompts
### Differential abundance and state
> "Run diffcyt DA with edgeR on my clustered SCE comparing treatment to control, and remind me how many biological replicates I need for this to be valid."
> "Run diffcyt DS with limma on the state markers within each cluster and show me which marker/cluster combinations change."

### Design
> "Set up a paired DA using a GLMM with patient as a random effect for my pre/post-treatment samples."
> "I have two batches confounded-ish with condition - put batch in the design matrix and explain when this can and can't be rescued."

### Compositionality
> "My regulatory-T cluster expanded a lot and now everything else looks depleted - re-test with a compositional method before I believe the depletion."

## What the Agent Will Do
1. Confirm biological replicates exist (>= 2-3 per group) and the sample is the unit.
2. Build a design + contrast (or a mixed-model formula for paired designs) from sample metadata.
3. Run diffcyt DA (edgeR/voom/GLMM) and/or DS (limma/LMM) via the `diffcyt()` wrapper on the SCE.
4. Apply FDR across clusters (and clusters x markers for DS).
5. Flag compositional artifacts and recommend simplex-aware re-testing when one population dominates.

## Tips
- The cell is not the unit - never run a per-cell test for a between-group claim; aggregate to per-sample summaries.
- Need >= 2-3 biological replicates per group; one sample per condition has no valid test.
- DA tests cluster frequency (type markers); DS tests within-cluster expression (state markers).
- Cluster proportions are compositional - a real expansion forces apparent depletion elsewhere; re-check dominant shifts with sccomp/scCODA/DCATS.
- Model batch in the design matrix; don't normalize it out then test naively. Fully confounded batch can't be rescued.
- Use `diffcyt-DA-GLMM` / `diffcyt-DS-LMM` with a random effect for paired/repeated-measures designs.
- diffcyt reuses edgeR/limma, so output semantics (logFC, p_adj) match differential-expression.

## Related Skills
- clustering-phenotyping - Cluster (on type markers) before testing
- gating-analysis - Compare manually gated population frequencies
- differential-expression/de-results - Shared edgeR/limma output semantics
- differential-expression/edger-basics - The count-model engine diffcyt reuses
- experimental-design/multiple-testing - FDR across clusters and clusters x markers
- experimental-design/batch-design - Model batch in the design, don't clean it out
