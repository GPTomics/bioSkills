# Quality Metrics - Usage Guide

## Overview

Gates IMC/MIBI data quality across pixel, channel, image, slide, and batch levels before analysis. The core point: QC is multi-level and every metric is blind at some level -- per-image SNR cannot see the failures that actually kill experiments (a dead antibody, an unbalanced batch, cells clustering by slide rather than phenotype). Counts are Poisson, dim is not failed, and IMC has no EQ-bead drift normalizer, so the discipline is to drop bad channels/ROIs/slides before analysis, not normalize after.

## Prerequisites

```bash
# Python
pip install numpy scikit-learn scanpy anndata matplotlib

# R / Bioconductor (spillover QC)
# BiocManager::install(c('CATALYST', 'spillR'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Compute cell-level SNR and tell me which channels failed versus which are just dim"
- "Read my spillover matrix and flag unacceptable crosstalk"
- "Check whether my cells cluster by patient instead of cell type"
- "Decide which ROIs to drop for striping or tissue loss"
- "Tell me if a marker is statistically indistinguishable from an empty channel"

## Example Prompts

### Channel QC
> "Fit a two-component mixture to each marker's per-cell counts and tell me which antibodies fail to separate positive from negative. Compare each to my 80ArAr channel."

### Spillover
> "My collagen channel is very bright. Which channels are at risk from its oxide and abundance-sensitivity spillover, and is the panel acceptable given my co-expression structure?"

### Drift and batch
> "I acquired over three weeks. How do I detect sensitivity drift without EQ beads, and how do I check for sample-of-origin clustering?"

### Gating decisions
> "Walk my ROIs and decide which to drop for striping, folding, or DNA dropout, and which regions to mask."

## What the Agent Will Do

1. Compute cell-level SNR via a two-component Gaussian mixture and compare each channel to a known-empty channel.
2. Read the spillover matrix by mass signature (M+-1 abundance, M+16 oxide, named-mass impurity), judging acceptability against co-expression.
3. Detect image-level artifacts (striping, folding, DNA/Ir dropout, saturation) and decide drop vs mask.
4. Monitor drift via the Ir/background channel across the slide and use an anchor reference sample for batch correction.
5. Embed cells and check for sample-of-origin clustering before any analysis; gate failures rather than normalizing them away.

## Tips

- A 1-2 count difference can be real biology -- do not borrow fluorescence SNR intuition.
- Dim is not failed; failure is inseparability of positive and negative populations, judged against an empty channel.
- High pixel correlation is not spillover; read the single-stain matrix and diagnose by mass signature.
- Compensation cannot rescue a saturated channel or separate real co-expression from leak; fix it in panel design.
- IMC has no in-line EQ-bead drift normalizer; use pre-tune criteria, watch the Ir channel, and anchor every batch with a reference sample.
- Always state thresholds in dual counts; many "standard" cutoffs are conventions, and some (in-line drift norm) do not exist.

## Related Skills

- data-preprocessing - hot-pixel removal, denoising, and NNLS spillover compensation
- cell-segmentation - segmentation QC and the impossible-co-expression monitor
- phenotyping - failed channels and batch corrupt cell-type calls
- differential-analysis - batch as a covariate when comparing across conditions
- flow-cytometry/cytometry-qc - suspension bead normalization and channel QC background
