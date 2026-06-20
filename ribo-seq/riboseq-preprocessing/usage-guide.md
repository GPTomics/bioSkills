# Ribo-seq Preprocessing - Usage Guide

## Overview

Preprocess ribosome profiling data from raw FASTQ: UMI extraction, adapter trimming, rRNA/tRNA contaminant depletion, footprint-aware alignment, deduplication (only when UMIs allow it), and read-length QC. The hard decisions are how the sample was harvested, which nuclease was used, whether to deduplicate, and which aligner fits the genome.

## Prerequisites

```bash
conda install -c bioconda cutadapt umi_tools sortmerna bowtie2 star samtools
```

## Quick Start

Tell your AI agent what you want to do:
- "Preprocess my Ribo-seq FASTQ files"
- "Extract UMIs and deduplicate my ribosome profiling reads"
- "Remove rRNA contamination from my footprints"
- "Align Ribo-seq reads with end-to-end settings"

## Example Prompts

### UMIs and Deduplication

> "My library has UMIs - extract them and deduplicate after alignment"

> "Should I deduplicate my Ribo-seq data? It has no UMIs"

> "Why are my high-occupancy codons flattened after marking duplicates?"

### Trimming and Contaminant Removal

> "Trim the 3' linker and discard reads with no adapter"

> "Remove rRNA contamination and report the percentage removed"

> "Build a combined rRNA plus tRNA index for depletion"

### Alignment and QC

> "Align my footprints with STAR using Ribo-seq-appropriate flags"

> "My data is from bacteria with MNase - how should I align it?"

> "Plot the read-length distribution and check for the 28-30 nt peak"

## What the Agent Will Do

1. Extract UMIs to the read name (if present) before any trimming
2. Trim the 3' linker with a permissive length floor and discard untrimmed reads
3. Deplete rRNA/tRNA contaminants before alignment
4. Align with STAR end-to-end (no soft-clipping) and transcriptome projection
5. Deduplicate on the BAM only when UMIs are present
6. QC the read-length distribution, contaminant fraction, and mapping rate

## Tips

- **Harvest method matters** - CHX pre-treatment distorts downstream dwell-time work; record drug and freezing
- **UMIs decide dedup** - with UMIs deduplicate; without UMIs do not (same position and length is mostly real)
- **rRNA is 50-90%** of raw reads - effective mRNA depth is a small fraction of the total
- **End-to-end alignment** - soft-clipping corrupts P-site offsets; the most important STAR change
- **Do not over-narrow size selection early** - plot the length distribution first, keep the 20-22 nt population
- **RNase I vs MNase** - bacteria need MNase and 3'-end anchoring; eukaryote-tuned settings break on bacterial data

## Related Skills

- ribosome-periodicity - Validate 3-nt periodicity and calibrate P-site offsets
- orf-detection - Detect translated ORFs from aligned footprints
- translation-efficiency - Requires matched, consistently processed RNA-seq
- read-qc/quality-reports - General read quality control
- read-alignment/star-alignment - General STAR alignment background
