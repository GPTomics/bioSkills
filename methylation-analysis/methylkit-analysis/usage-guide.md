# methylKit Analysis - Usage Guide

## Overview
methylKit is the R/Bioconductor object model for short-read bisulfite (WGBS/RRBS) differential methylation. This skill owns the import-to-results spine: reading Bismark coverage or cytosine-report files into a methylRawList, then filtering by coverage, normalizing, uniting (with optional destranding), and testing for differential methylation at single CpGs (DMCs) or fixed tiles (DMRs). It also covers PCA/correlation/clustering QC and batch handling. The point that sets your results: filtering, normalization, destranding, and the overdispersion model are decisions that set the false-positive rate and which sites survive - the defaults offer no protection.

## Prerequisites
```r
if (!require('BiocManager', quietly=TRUE)) install.packages('BiocManager')
BiocManager::install(c('methylKit', 'GenomicRanges', 'genomation', 'annotatr'))
```

Conceptual prerequisites:
- Bismark coverage (`.cov`/`.cov.gz`) or cytosine report (`*.CpG_report.txt`) input - produced upstream by methylation-calling. Destranding requires the cytosine report (the `.cov` format carries no strand).
- Biological replicates per group. Overdispersion correction and any valid inference need replicates; pooling them away invalidates the test.
- The genome build of the coverage files must match any annotation package used downstream; the `assembly=` string is metadata only.
- Large WGBS objects can exhaust memory - use `dbtype='tabix'`, `save.db=TRUE` for disk-backed methylRawDB objects.

## Quick Start
Tell your AI agent what you want to do:
- "Load my Bismark coverage files into methylKit and run QC"
- "Filter low-coverage CpGs, normalize, and unite my samples"
- "Find differentially methylated CpGs with overdispersion correction"
- "Run a fast tile-based DMR screen with methylKit"
- "Check whether my samples cluster by treatment or by batch"

## Example Prompts

### Data Import
> "Read my four Bismark coverage files into methylKit as two controls and two treated samples, then show coverage and methylation stats."

> "Import my Bismark cytosine reports so I can destrand CpGs properly."

### Filtering and Normalization
> "Filter CpGs below 10x and above the 99.9th coverage percentile, then median-normalize coverage across samples before uniting."

### Quality Control
> "Generate PCA, sample correlation, and hierarchical clustering for my united methylation object and tell me whether samples separate by treatment."

### Differential Analysis
> "Test individual CpGs for differential methylation with MN overdispersion correction and BH adjustment, then report DMCs with at least a 25% difference at q < 0.01, split into hyper and hypo."

> "Run a tile-based DMR screen at 1 kb windows requiring at least 3 covered CpGs per tile, and flag that the per-tile q is a screen, not a selection-corrected region FDR."

### Batch and Design
> "My PC1 tracks sequencing batch - associate components with the batch covariate and remove the batch component before testing."

## What the Agent Will Do
1. Build sample metadata (file paths, sample IDs, integer treatment vector) and import with methRead, choosing the pipeline that matches the input format.
2. Generate per-sample coverage and methylation QC stats.
3. Filter by coverage (low tail and high tail) and normalize coverage between samples - before uniting or tiling.
4. Unite samples into a per-base methylBase object, destranding only for CpG with strand information present.
5. QC the united object with correlation, PCA, and clustering to check for batch confounding.
6. Run calculateDiffMeth with overdispersion='MN' (which uses the F-test) and an explicit adjust method.
7. Filter to significant DMCs with getMethylDiff by effect size and q-value, separating hyper and hypo.
8. Optionally tile (filtered/normalized) into windows with cov.bases >= 3 for a fast region screen, noting the q is not selection-corrected.
9. Export results to data frame or BED and hand off to annotation/enrichment.

## Output Interpretation

| Column | Description |
|--------|-------------|
| chr, start, end | Genomic position of the CpG or tile |
| meth.diff | Methylation difference (%); positive = hyper in the higher-treatment group |
| pvalue | Raw p-value (F-test under MN, otherwise the chosen test) |
| qvalue | Adjusted p-value (SLIM by default; BH only if set) |

## Tips
- Use `pipeline='bismarkCoverage'` for `.cov` files and `pipeline='bismarkCytosineReport'` for CX/CpG reports; only the report enables proper destranding.
- Set `overdispersion='MN'` whenever replicates exist; the default 'none' over-calls. MN forces the F-test, so do not also pass `test='Chisq'` and expect a chi-square.
- The default `adjust` is SLIM, not BH. Set `adjust='BH'` if the q-values must be comparable to DSS, limma, or another tool.
- Filter and normalize BEFORE `unite` and BEFORE `tileMethylCounts`; tiling the raw object propagates artifacts.
- Raise `tileMethylCounts(cov.bases=)` to at least 3 so a window needs real CpG support; the default 0 admits single-CpG "regions."
- The 25% difference and q < 0.01 thresholds are conventions, not derived values - justify them by feature, coverage, and cell purity, and report them.
- `pool()` destroys biological replication; use it only for no-replicate exploratory visualization, never for the reported test.
- methylKit tiles are a fast screen; for selection-aware region FDR use dmr-detection (dmrseq).

## Related Skills

- methylation-calling - Produces the coverage/cytosine reports read here
- differential-cpg-testing - Per-site statistical model choice (count vs continuous)
- dmr-detection - Selection-aware region callers (dmrseq/DSS) beyond methylKit tiles
- pathway-analysis/go-enrichment - Functional annotation of differentially methylated genes
- long-read-sequencing/nanopore-methylation - Long-read MM/ML calling; pipe counts into this object model
- workflows/methylation-pipeline - End-to-end bisulfite pipeline
