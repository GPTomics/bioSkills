# ATAC-seq Quality Control - Usage Guide

## Overview

Compute the seven canonical ATAC-seq QC metrics (depth, alignment rate, mitochondrial fraction, library complexity, fragment-size periodicity, TSS enrichment, FRiP) and grade each against ENCODE 4 acceptance thresholds. Diagnose specific failure modes (over-transposition, library bottlenecking, chromatin degradation, batch effects) and decide which replicates to drop before peak calling.

## Prerequisites

```bash
conda install -c bioconda samtools picard deeptools multiqc bedtools
pip install pysam pyBigWig numpy pandas matplotlib
```

```r
BiocManager::install(c('ATACseqQC', 'TxDb.Hsapiens.UCSC.hg38.knownGene', 'BSgenome.Hsapiens.UCSC.hg38'))
```

ENCODE blacklist (Amemiya 2019), GENCODE TSS BED (matched to genome build) required for TSS enrichment.

## Quick Start

Tell your AI agent what you want to do:
- "Score every ENCODE 4 ATAC QC metric and flag PASS/WARN/FAIL"
- "Compute TSS enrichment using the ENCODE pyTSSe formula, not ATACseqQC's"
- "Diagnose why fragment-size distribution is flat -- over-transposition vs degraded chromatin"
- "Calculate NRF, PBC1, and PBC2 from the raw mapped BAM (pre-dedup)"
- "Run plotFingerprint and compare fingerprint QC metrics across replicates (the plain JS distance column needs a --JSDsample reference)"

## Example Prompts

### ENCODE 4 Compliance Audit
> "Run a full ENCODE 4 ATAC-seq QC audit on three replicates: nuclear read depth, alignment rate, mt fraction, NRF/PBC1/PBC2, TSS enrichment (ENCODE method), FRiP using MACS peaks, and fragment-size periodicity. Flag any metric that fails ENCODE thresholds and recommend drop/keep per replicate."

### TSS Enrichment Reconciliation
> "My ATACseqQC TSSEscore is 21 but my collaborator says ENCODE wants >= 7. Compute the ENCODE pyTSSe-style score so we can compare apples to apples."

### Diagnosing Flat Fragment Distribution
> "Insert-size distribution shows one broad peak with no nucleosome periodicity. Diagnose the likely cause (over-transposition, degraded chromatin, MNase contamination) and tell me whether the library is salvageable."

### Library Complexity Investigation
> "NRF is 0.65 and PBC2 is 0.8. Confirm bottlenecking vs legitimate ATAC duplicates at hyperaccessible sites by checking pile-up coverage at top peaks."

### Replicate Reconciliation
> "Pearson correlation between rep2 and rep3 is 0.78 (others are 0.94). Run PCA on peak counts and check whether rep2 and rep3 cluster together by batch -- if so, plan to add batch as a covariate downstream."

### Cross-Study Comparison
> "Compare my TSS enrichment to a published study; we both used hg38 but different TSS sources. Recompute on the same GENCODE v29 protein-coding TSS file for parity."

## What the Agent Will Do

1. Run `samtools flagstat` and `samtools idxstats` for alignment summary and chrM count
2. Compute NRF / PBC1 / PBC2 on the pre-dedup mapped BAM (post-MAPQ-filter)
3. Generate fragment-size distribution with Picard CollectInsertSizeMetrics; classify periodicity pattern
4. Compute TSS enrichment (ENCODE pyTSSe convention by default; ATACseqQC if requested)
5. Calculate FRiP from MACS peaks
6. Run plotFingerprint (synthetic JS distance; the plain JS distance requires a --JSDsample reference); cross-replicate Spearman correlation via `multiBamSummary`
7. Aggregate everything via MultiQC
8. Grade each metric PASS/WARN/FAIL and write a per-sample report card

## ENCODE 4 Threshold Quick Reference

| Metric | Ideal | Acceptable | Reject |
|--------|-------|------------|--------|
| Nuclear reads (post-dedup, no chrM) | >= 50M | 25-50M | < 25M |
| Alignment rate | >= 95% | 80-95% | < 80% |
| Mt fraction (Omni-ATAC) | < 5% | < 20% | > 50% |
| NRF | >= 0.9 | 0.7-0.9 | < 0.7 |
| PBC1 | >= 0.9 | 0.7-0.9 | < 0.7 |
| PBC2 | >= 3.0 | 1.0-3.0 | < 1.0 |
| TSS enrichment (ENCODE) | >= 7 | 5-7 | < 5 |
| FRiP | >= 0.3 | 0.2-0.3 | < 0.2 |

## Fragment-Size Pattern Quick Reference

| Pattern | Diagnosis |
|---------|----------|
| NFR + mono + di + tri visible | Excellent |
| NFR + mono only, di faint | Acceptable (Omni-ATAC) |
| Single broad peak, no NFR | Over-transposition; reject |
| Mono >> NFR | Under-transposition; caution |
| Sharp 147 bp spike, no flanks | MNase-like; reject |
| 10.4 bp helical sub-peaks visible | High-quality structural resolution |

## Tips

- TSS enrichment formula MUST be specified: ENCODE pyTSSe and ATACseqQC TSSEscore differ by 2-3x because of window size choices. State which implementation was used.
- Compute NRF/PBC on the raw mapped BAM (after MAPQ filter, before dedup). Computing post-dedup gives NRF = 1.0 trivially.
- ENCODE mt-fraction thresholds assume Omni-ATAC. Standard ATAC-seq routinely shows 30-50% chrM and is still usable after chrM-stripping; mt > 50% is the practical reject line.
- NRF below 0.9 is not by itself fatal: ATAC has *legitimate* duplicates at hyperaccessible sites. Confirm bottlenecking with PBC1 < 0.7 + PBC2 < 1.0 + visual pile-ups.
- FRiP is peak-set-dependent. Pin the peak caller version and q-value cutoff before reporting FRiP across studies.
- Use Spearman correlation (not Pearson) for cross-replicate comparison; ATAC signal is heavy-tailed and Pearson is dominated by top peaks.
- The 10.4 bp helical periodicity in fragment-size distribution is a positive QC indicator, not a requirement. Its absence does not flag the library.
- For non-model organisms, ENCODE thresholds do not apply directly. Use cohort percentile rank (e.g., drop the bottom 10% of TSS enrichment within the cohort).
- multiBamSummary `--skipZeros` is mandatory for ATAC; otherwise empty regions dominate the correlation.

## Related Skills

- atac-seq/atac-peak-calling - QC informs which replicates to feed peak calling
- atac-seq/nucleosome-positioning - Detailed fragment-size analysis and V-plots
- atac-seq/single-cell-atac - per-cell QC with different thresholds
- read-qc/quality-reports - Upstream FASTQ QC
- alignment-files/bam-statistics - samtools flagstat / idxstats wrappers
- alignment-files/duplicate-handling - Pre-dedup before NRF/PBC computation
