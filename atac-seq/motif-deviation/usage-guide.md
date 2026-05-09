# Motif Deviation (chromVAR) - Usage Guide

## Overview

Compute per-sample (or per-cell) TF motif accessibility z-scores using chromVAR. Identifies which motifs covary with sample/cell state, controlling for GC content and overall accessibility through matched background peak sets. Complements footprinting (which scores individual binding sites) with genome-wide motif-class summary statistics.

## Prerequisites

```r
BiocManager::install(c('chromVAR', 'motifmatchr', 'JASPAR2024', 'TFBSTools',
                       'BSgenome.Hsapiens.UCSC.hg38', 'SummarizedExperiment',
                       'limma', 'Signac', 'Seurat'))

# ArchR alternative (if using ArchR ecosystem)
remotes::install_github('GreenleafLab/ArchR', ref='master', repos=BiocManager::repositories())
```

Inputs:
- Peak count matrix (rows = peaks, cols = samples) OR Seurat/ArchR object
- Genomic ranges of peaks (BED or GRanges)
- Motif PFM database (JASPAR 2024 default; HOCOMOCO or CIS-BP for special cases)

## Quick Start

Tell your AI agent what you want to do:
- "Run chromVAR on a 12-sample bulk ATAC peak count matrix and report top variable motifs"
- "Run chromVAR per-cell on Signac scATAC object via RunChromVAR"
- "Compare motif z-scores between treatment groups using limma on the deviation matrix"
- "Use ArchR addDeviationsMatrix on PeakMatrix and rank cluster-marker motifs"
- "Verify chromVAR has enough peaks (>5000) and reads (>1500/sample) for stable inference"

## Example Prompts

### Bulk Analysis
> "Run chromVAR on this peak count matrix with JASPAR 2024 vertebrate CORE motifs. Add GC bias correction, filter samples below 1500 reads-in-peaks and FRiP < 0.15, then report top 20 most variable motifs and their z-score matrix."

### Differential Motif Activity
> "Run limma on the chromVAR z-score matrix to identify motifs differing between control and treated. Report adj.P.Val < 0.05 and abs(logFC) > 0.5; do not use raw counts."

### Single-Cell Signac
> "Add motifs to my Seurat object after the peakset is finalized, run RunChromVAR, then call FindAllMarkers on the chromvar assay with `mean.fxn=rowMeans` and `fc.name='avg_diff'` to find cluster-marker motifs."

### Single-Cell ArchR
> "In my ArchR project, add reproducible peakset, addPeakMatrix, addMotifAnnotations with `motifSet='cisbp'`, addBgdPeaks, addDeviationsMatrix. Then getMarkerFeatures on MotifMatrix grouped by Clusters."

### Validation
> "My chromVAR variability scores are all >5; verify peak count is at full ATAC scale (50k+) and per-sample reads-in-peaks > 1500. If too sparse, aggregate cells before running."

### Time-Course
> "Fit a spline regression on chromVAR z-scores across the time course; identify motifs with non-monotonic trajectories indicative of waves of TF activity."

## What the Agent Will Do

1. Verify peakset (>5000 peaks) and per-sample depth (>1500 reads-in-peaks for bulk)
2. Build SummarizedExperiment from counts + peak ranges
3. Add GC bias annotation (chromVAR requires genomic sequence)
4. Filter samples and peaks per chromVAR thresholds
5. Match motifs (JASPAR 2024 default; CIS-BP for non-model)
6. Sample matched background peaks (50 iterations; default is fine)
7. Compute deviation matrix (motif x sample)
8. Compute variability (per-motif variance across samples)
9. Optionally: differential testing via limma; per-cluster markers via FindMarkers / getMarkerFeatures
10. Generate visualization: variability plot, heatmap of top motifs, PCA

## Tool Decision Quick Reference

| Setting | Tool |
|---------|------|
| Bulk, 6+ samples | chromVAR + limma |
| Bulk, 3-5 samples | chromVAR; variability ranking only |
| scATAC, Seurat | Signac::RunChromVAR -> FindAllMarkers |
| scATAC, ArchR | addDeviationsMatrix -> getMarkerFeatures |
| Multiomics scATAC + scRNA | chromVAR + paired DE; or SCENIC+ |
| Plant / non-model | chromVAR with custom CIS-BP motifs and BSgenome |

## Variability Interpretation Quick Reference

| Variability | Interpretation |
|-------------|----------------|
| < 1 | Motif activity ~constant |
| 1-2 | Modest variation |
| 2-5 | Strong variation; biologically interesting |
| > 5 | Major driver; flagship hit |

## chromVAR vs Footprinting Decision

| Question | Use |
|----------|-----|
| Which TFs distinguish my conditions/clusters? | chromVAR |
| Is THIS specific motif site bound? | Footprinting (TOBIAS/HINT) |
| Per-cell TF activity for trajectory analysis | chromVAR (Signac/ArchR) |
| Motif activity correlated with paired RNA | chromVAR + co-expression OR SCENIC+ |

## Tips

- chromVAR requires biological variation across samples to be informative. With only one condition replicated, use footprinting or differential accessibility instead.
- Minimum thresholds: 5000+ peaks, 1500+ reads-in-peaks per sample (bulk), 500+ cells per cluster (single-cell). Below these, z-scores are unreliable.
- Default `niterations=50` for `getBackgroundPeaks` is well-calibrated; reducing it adds noise, increasing it doesn't help. Don't tune unless benchmarking.
- chromVAR z-scores are NOT comparable across studies that use different peak sets or different motif databases. Recompute for cross-study comparison.
- Use BH-adjusted p-values from limma (`adj.P.Val`) for differential motif testing on z-scores; raw t-tests are mis-calibrated.
- For single-cell, run `RunChromVAR()` AFTER finalizing the peakset; otherwise motif annotations become stale.
- ArchR's `cisbp` motif set has ~5000 motifs (broader); JASPAR has ~1900 (more curated). Choose based on whether breadth or precision is needed.
- chromVAR is fundamentally a *summary statistic*. For per-site analysis, use footprinting; for differential peak counts, use DiffBind.
- The variability score does not require condition labels and is a good unsupervised discovery tool. Z-scores require labels for differential.
- When global accessibility differs >5x between cell types in the dataset, run chromVAR per-cluster separately rather than jointly.

## Related Skills

- atac-seq/footprinting - Per-site TF binding analysis
- atac-seq/differential-accessibility - Peak-level differential
- atac-seq/single-cell-atac - sc workflow integration
- atac-seq/co-accessibility - Cis-regulatory follow-up
- gene-regulatory-networks/scenic-regulons - Downstream TF -> target networks
- chip-seq/motif-analysis - Alternative motif enrichment
- single-cell/clustering - Cluster definition for per-cluster analysis
