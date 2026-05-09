# ATAC-seq Peak Calling - Usage Guide

## Overview

Call accessible chromatin regions from ATAC-seq BAM files. Covers MACS2/MACS3 (ENCODE default), Genrich (joint replicates), MACS3 hmmratac (HMM-based), and downstream IDR for reproducible peak sets. Adapts ENCODE 4 ATAC-seq pipeline conventions for shift-extend modeling, pseudoreplicate consistency, and blacklist/greylist filtering.

## Prerequisites

```bash
conda install -c bioconda macs3 macs2 genrich samtools bedtools idr
# HMMRATAC is now bundled in MACS3 as `macs3 hmmratac` (since MACS 3.0.0)

# Blacklist files
wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
gunzip hg38-blacklist.v2.bed.gz
```

Inputs assumed deduplicated and chrM-removed. See `alignment-files/duplicate-handling` and `read-alignment/bowtie2-alignment` for the upstream alignment + filtering steps.

## Quick Start

Tell your AI agent what you want to do:
- "Call ATAC-seq peaks following the ENCODE 4 pipeline with IDR across replicates"
- "Use Genrich joint mode to call peaks across two replicates and exclude chrM"
- "Run HMMRATAC on a deep ATAC library to get NFR + flanking nucleosome calls"
- "Call peaks on nucleosome-free fragments only for footprinting input"
- "Reproduce peaks from a published dataset using the same effective genome size as deepTools"

## Example Prompts

### ENCODE 4 Pipeline
> "Run the ENCODE 4 ATAC-seq peak-calling pipeline: per-replicate MACS2 with `-p 0.01`, then pseudoreplicate split + IDR with `--idr-threshold 0.10` for self-consistency, and IDR `<= 0.05` for true replicates. Apply ENCODE Nself <= 2 rule to flag failing libraries."

### Joint Replicate Calling with Genrich
> "Call ATAC peaks jointly across three replicates with Genrich, removing chrM (`-e chrM`), filtering ENCODE blacklist (`-E hg38-blacklist.v2.bed`), and removing PCR duplicates inside Genrich (`-r`)."

### HMM-based Calling
> "Run `macs3 hmmratac` on a 50M-read ATAC library to produce open / nucleosomal / background calls and verify fragment-size periodicity is strong enough first."

### NFR-only for Footprinting
> "Filter the BAM to fragments < 100 bp and call peaks with `--shift -37 --extsize 75` to produce a footprinting-ready peak set."

### Single-Sample (No Replicate)
> "Call peaks on a single ATAC sample with `-q 0.05` (tighter than ENCODE) since IDR is not applicable, then filter against the ENCODE blacklist."

### Effective Genome Size
> "Look up the deepTools effectiveGenomeSize for hg38 at our 100 bp read length and use that as `-g` instead of the MACS shorthand."

## What the Agent Will Do

1. Verify upstream BAM is deduplicated, MAPQ-filtered, and chrM-stripped
2. Choose caller based on depth, replicate count, and downstream use (footprinting vs differential vs domain-level)
3. Pick `-f BAM` (with `--shift/--extsize`) for ENCODE-style or `-f BAMPE` (no shift) for fragment-aware mode
4. Set effective genome size from deepTools table, not the legacy `-g hs/mm` shorthand
5. Per-replicate, pooled, and pseudoreplicate peak calls if running ENCODE pipeline
6. Run IDR on true replicates (threshold 0.05) and pseudoreplicates (threshold 0.10)
7. Apply ENCODE Nt/Nself <= 2 self-consistency rule
8. Filter ENCODE blacklist and optionally a sample-derived greylist
9. Produce narrowPeak + bigWig signal track outputs

## Caller Decision Quick Reference

| Situation | Caller |
|-----------|--------|
| 2-3 replicates, depth >= 25M | MACS2 ENCODE pipeline + IDR |
| 1 sample, no replicates | MACS3 callpeak `-q 0.05` (no IDR) |
| Want NFR + flanking nucleosomes | MACS3 hmmratac |
| Joint multi-rep with built-in chrM/blacklist | Genrich `-j -e chrM -E blacklist` |
| Broad super-enhancer regions | MACS3 `--broad --broad-cutoff 0.1` |
| FFPE / degraded chromatin | MACS3 callpeak (avoid HMM) |
| scATAC pseudobulk | MACS3 per cluster, see single-cell-atac |
| Plant / non-model | MACS3 with empirically computed `-g` |

## Effective Genome Size Quick Reference

| Genome | 50 bp reads | 75 bp reads | 100 bp reads |
|--------|-------------|-------------|--------------|
| hg38 | 2.913e9 | 2.747e9 | 2.701e9 |
| hg19 | 2.864e9 | 2.770e9 | 2.701e9 |
| mm10 | 2.652e9 | 2.467e9 | 2.407e9 |
| mm39 | 2.654e9 | 2.494e9 | 2.494e9 |

Source: deepTools `effectiveGenomeSize` documentation. Values change by read length because shorter reads have more multimappers excluded.

## ENCODE Self-Consistency Rule

Library passes if `max(Nt, Nself) / min(Nt, Nself) <= 2` where:
- **Nt** = peaks passing IDR <= 0.05 on true replicate pair
- **Nself** = peaks passing IDR <= 0.10 on pseudoreplicate pair (single rep split in half)

Both ratios > 2 means the library is rejected per ENCODE 4 standards.

## Tips

- ENCODE pipeline uses MACS2 with `-p 0.01` (loose, for IDR), not `-q 0.05`. Tighten only if NOT running IDR.
- `-f BAMPE` silently ignores `--shift/--extsize`. To use shift modeling, use `-f BAM` and treat ends as single reads.
- `--keep-dup all` is mandatory for ATAC: Tn5 generates legitimate duplicate cuts at hyperaccessible sites.
- Always remove chrM before peak calling (`samtools view -h sample.bam | grep -v "chrM"`); chrM accumulates ATAC reads at >50% of total in poor preps.
- For broad accessibility regions (super-enhancers, MYOD1 regulons), narrow mode fragments them; use `--broad --broad-cutoff 0.1`. But do not run IDR on broad peaks.
- Re-center peaks on summits +/- 250 bp (Corces 2017 Omni-ATAC convention) for differential analysis to avoid width-driven count differences.
- IDR ranking column matters: use `--rank p.value` (column 8 of narrowPeak); `--rank signal.value` is unreliable when MACS pileup scaling differs.
- HMMRATAC needs >= 30M deduplicated nuclear reads and clear fragment-size periodicity; verify with the atac-qc skill first.
- Genrich's joint mode does NOT need pre-deduplication when `-r` is used; doing both is double-removing.
- The ENCODE blacklist (Amemiya 2019 v2) is mandatory; greylist is optional and most useful for matched same-batch libraries.

## Related Skills

- atac-seq/atac-qc - QC checks before peak calling (TSS enrichment, fragment periodicity)
- atac-seq/consensus-peakset - Build differential-ready fixed-width peaks
- atac-seq/single-cell-atac - Pseudobulk peak calling per cluster
- atac-seq/differential-accessibility - Differential testing on the peak set
- atac-seq/footprinting - Use NFR-only peak calls as footprinting input
- read-alignment/bowtie2-alignment - Upstream ATAC alignment
- alignment-files/duplicate-handling - Pre-call BAM filtering
- chip-seq/peak-calling - Compare to ChIP-seq peak calling (uses input control)
- genome-intervals/bed-file-basics - narrowPeak/BED file operations
