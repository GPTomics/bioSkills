# Consensus Peakset Construction - Usage Guide

## Overview

Build a single non-redundant peakset that all samples can be counted against. Covers the modern ATAC standards: Corces 2018 iterative overlap removal (501 bp fixed-width), DiffBind summit re-centering, ENCODE IDR-based consistency, and per-condition union strategies. Width and overlap-rule choices propagate to every downstream analysis (differential, motif enrichment, ML feature engineering).

## Prerequisites

```bash
conda install -c bioconda bedtools samtools subread bedops
pip install pybedtools
```

```r
BiocManager::install(c('DiffBind', 'GenomicRanges'))
```

Inputs: per-replicate peak files (preferably narrowPeak with summit info in column 10), genome chrom sizes, ENCODE blacklist (v2).

## Quick Start

Tell your AI agent what you want to do:
- "Build a Corces 2018 iterative-overlap consensus peakset (501 bp fixed-width) from per-replicate narrowPeak files"
- "Use DiffBind dba.count summits=250 for fixed-width counting"
- "Per-condition consensus then union to preserve condition-specific peaks"
- "ENCODE-style IDR-based consensus across rep pairs at threshold 0.05"
- "Filter consensus against hg38 blacklist v2 and convert to SAF for featureCounts"

## Example Prompts

### Modern Iterative Overlap (Corces 2018)
> "Pool all per-replicate narrowPeak files; re-center each peak on its summit (column 10 offset); extend +/- 250 bp; sort by signalValue descending; greedily keep non-overlapping peaks. Filter ENCODE blacklist v2 and report final peak count."

### DiffBind Summit Re-centering
> "Use DiffBind dba.count with `summits=250` (501 bp fixed-width) and `minOverlap=2` (peak in >= 2 reps) on this 6-sample dataset."

### Per-Condition Strategy
> "Build a consensus per condition first (peak in >= 2/3 reps within condition), then union across conditions. This preserves condition-specific peaks that majority rule would lose."

### ENCODE IDR
> "Run IDR per rep pair at threshold 0.05 (true reps); union the IDR-passed peaks; filter blacklist; convert to SAF for featureCounts."

### featureCounts Matrix Generation
> "Convert the consensus BED to SAF format with chr_start_end as GeneID; run featureCounts with `-p --countReadPairs` on all sample BAMs to produce the count matrix for DESeq2."

### Cross-Study Consensus
> "I have peaks from a published study (hg19) and my own (hg38). Lift over the published peaks to hg38, then merge with mine using iterative overlap to enable cross-study comparison."

## What the Agent Will Do

1. Verify per-replicate peaks have summit info (narrowPeak column 10) for re-centering
2. Choose strategy based on goal: iterative overlap for ML/cross-study, DiffBind for standard DA, per-condition union for strong biology shifts
3. Re-center on summits +/- 250 bp (standard) or other half-width
4. Apply overlap rule (greedy non-overlap, IDR-pass, majority rule)
5. Filter ENCODE blacklist v2
6. Convert to SAF for featureCounts OR keep as BED for DiffBind dba.count
7. Document the strategy explicitly in the methods

## Strategy Decision Quick Reference

| Goal | Strategy |
|------|----------|
| ML feature engineering / cross-study | Iterative overlap (Corces 2018), 501 bp fixed |
| Standard DiffBind DA workflow | DiffBind summits=250, minOverlap=2 |
| ENCODE-compliant differential | IDR-pass per pair, true-rep threshold 0.05 |
| Strong condition-specific biology | Per-condition consensus then union |
| scATAC pseudobulk per cluster | MACS3 per cluster + iterative overlap across clusters |
| Quick exploratory | bedtools merge (acknowledge width bias) |

## Width Choice Quick Reference

| Half-width | Total | Use |
|-----------|-------|-----|
| 100 bp | 201 bp | Fine-resolution; sub-peak structure preserved |
| 150 bp | 301 bp | Mid-resolution |
| 250 bp | 501 bp | Corces 2018 / DiffBind default |
| 500 bp | 1001 bp | Coarse; broad regulatory regions |

501 bp is the modern de facto standard; smaller widths preserve resolution at the cost of more peaks; larger widths consolidate but lose specificity.

## Tips

- Always document the strategy in methods. Strategy choice can shift downstream peak count and FDR by 30-50%.
- Iterative overlap is the modern standard for ATAC ML feature engineering. DiffBind summits=250 is the default for standard differential workflows. Both are valid; they answer the same biology slightly differently.
- Variable-width consensus drives width-confounded differential. Always re-center to fixed width.
- Per-condition consensus then union preserves condition-specific peaks; global majority rule discards them.
- Filter ENCODE blacklist v2 AFTER consensus construction, BEFORE counting. Filtering before consensus loses some IDR-passed peaks.
- Summit position is narrowPeak column 10 (offset from start); use `start + summit_offset` as the absolute summit. broadPeak files don't have this; midpoint approximation is acceptable but biased.
- IDR thresholds differ for true reps (0.05) vs pseudoreps (0.10). Mix carefully: pseudorep IDR is looser and inflates the consensus.
- featureCounts SAF requires unique GeneID; use chr_start_end naming to guarantee uniqueness.
- For paired-end ATAC, featureCounts must use `-p --countReadPairs` (counts pairs, not individual reads).
- For peak-to-gene assignment downstream, the consensus peakset is annotated once via ChIPseeker; results stay stable across reruns.

## Related Skills

- atac-seq/atac-peak-calling - Generate per-replicate peaks (input)
- atac-seq/differential-accessibility - Use consensus for DA testing
- atac-seq/atac-qc - Filter failing samples before consensus
- atac-seq/single-cell-atac - Per-cluster consensus for pseudobulk
- genome-intervals/bed-file-basics - bedtools operations
- genome-intervals/interval-arithmetic - merge / intersect / subtract
- chip-seq/peak-calling - Same consensus strategy for ChIP
