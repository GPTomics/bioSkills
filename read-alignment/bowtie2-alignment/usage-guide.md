# Bowtie2 Alignment - Usage Guide

## Overview

Bowtie2 is a fast, memory-efficient short-read aligner whose defining choice is end-to-end (the whole read must align) vs local (read ends may be soft-clipped) mode. It is the de-facto aligner for ChIP-seq, ATAC-seq, and CUT&RUN, and it is the engine Bismark wraps for bisulfite data. For peak assays the alignment itself matters less than the mode, the sensitivity preset, and the fragment-geometry flags (--no-mixed, --no-discordant, --dovetail, -X) that determine the fragment coordinates the peak caller consumes -- this skill emphasizes those.

## Prerequisites

```bash
conda install -c bioconda bowtie2 samtools
```

- A Bowtie2 index (built with `bowtie2-build`) or a prebuilt index; pass the basename to `-x`.
- For ChIP/ATAC, decide whether reads are trimmed (end-to-end) or carry adapter read-through (local) before choosing the mode.

## Quick Start

Tell your AI agent what you want to do:
- "Align my ChIP-seq reads with Bowtie2 and filter multimappers"
- "Build a Bowtie2 index for my reference genome"
- "Align ATAC-seq reads with the right fragment and soft-clipping settings"
- "Choose end-to-end vs local mode for my data and explain why"
- "Report up to 5 alignments per read for a repetitive region"

## Example Prompts

### ChIP-seq
> "Align my paired-end ChIP-seq reads with high sensitivity, keep only concordant proper pairs, and filter to MAPQ >= 30 before peak calling."

### ATAC-seq
> "Align ATAC-seq reads allowing fragments up to 2 kb and soft-clipping adapter read-through, and tell me which flags handle the short dovetailed fragments."

### Mode choice
> "My reads have adapter contamination at the 3' end and a low alignment rate in end-to-end mode -- should I trim or switch to local mode?"

### Multi-mapping
> "Align reads and report up to 5 alignments per read for a repetitive region, and explain why MAPQ is unreliable in -k mode."

## What the Agent Will Do

1. Build the index with `bowtie2-build` if needed, and pass the basename (not a filename) to `-x`.
2. Choose end-to-end vs local mode from the data (trimmed clean DNA vs adapter-contaminated/ATAC) and a sensitivity preset.
3. Set the fragment-geometry flags for the assay: `--no-mixed --no-discordant` for ChIP, plus `--local --dovetail -X 2000` for ATAC.
4. Filter to a Bowtie2-appropriate MAPQ (`-q 30`, not a BWA-style `-q 60`) and a flag mask (`-F 1804`) before peak calling.
5. Add read groups (`--rg-id`/`--rg`) and stream to a sorted, indexed BAM.
6. Route peak calling and the ATAC Tn5 cut-site shift to the chip-seq / atac-seq categories, and the QC gate to alignment-files/bam-statistics.

## When to Use Bowtie2 vs bwa-mem2

| Use case | Bowtie2 | bwa-mem2 |
|----------|---------|----------|
| ChIP-seq / CUT&RUN | preferred (ENCODE standard) | works |
| ATAC-seq | preferred (local + dovetail) | works |
| DNA variant calling (WGS/WES) | not the community default | preferred |
| RNA-seq | no (use STAR/HISAT2) | no (use STAR/HISAT2) |
| Soft-clipping adapter-contaminated ends | yes (`--local`) | yes (`-L` tuning) |

## Tips

- Choose the mode from the data: end-to-end for clean trimmed DNA, `--local` for adapter read-through or ATAC.
- Never copy a BWA `MAPQ >= 60` filter to Bowtie2 output -- the scale caps at 42 (end-to-end) / 44 (local); use `-q 30`.
- For ATAC, add `--dovetail` and `-X 2000` or real short-fragment pairs get flagged discordant and dropped.
- Pass the index basename to `-x`, not a `.bt2` file.
- Bowtie2 writes its alignment summary to stderr -- capture it (`2> align.log`) for QC and MultiQC aggregation.
- The Tn5 +4/-5 cut-site shift for ATAC is a downstream signal-track transform, not an alignment flag; route it to the atac-seq category.

## Related Skills

- bwa-alignment - DNA variant-calling alignment with bwa-mem2 (ALT/decoy-aware)
- star-alignment - RNA splice-aware alignment (when reads cross junctions)
- read-qc/fastp-workflow - Trim adapters before end-to-end alignment
- alignment-files/duplicate-handling - Mark/remove duplicates after alignment
- alignment-files/bam-statistics - flagstat/idxstats QC gate; the cross-tool MAPQ scale via sam-bam-basics
- chip-seq/peak-calling - Call peaks from ChIP/CUT&RUN BAMs
- atac-seq/atac-peak-calling - ATAC peak calling and the Tn5 cut-site shift
- methylation-analysis/bismark-alignment - Bisulfite alignment (wraps Bowtie2)
