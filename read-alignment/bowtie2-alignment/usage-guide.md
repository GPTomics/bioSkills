# Bowtie2 Alignment Usage Guide

## Overview

Bowtie2 is a fast and memory-efficient aligner for short reads. It supports both end-to-end and local alignment modes, making it versatile for various applications including ChIP-seq, ATAC-seq, and general short-read alignment.

## Installation

```bash
conda install -c bioconda bowtie2 samtools
```

## Quick Start

```bash
# 1. Build index
bowtie2-build reference.fa bt2_index

# 2. Align reads
bowtie2 -p 8 -x bt2_index -1 r1.fq.gz -2 r2.fq.gz | samtools sort -o aligned.bam -
samtools index aligned.bam
```

## When to Use Bowtie2 vs BWA

| Use Case | Bowtie2 | BWA-MEM2 |
|----------|---------|----------|
| ChIP-seq | Preferred | Works |
| ATAC-seq | Preferred | Works |
| RNA-seq | No (use STAR) | No |
| WGS | Works | Preferred |
| Local alignment | Yes | Limited |

## Alignment Mode Selection

### End-to-End (Default)
Aligns the entire read, penalizing unaligned bases:
- Best for: High-quality reads, variant calling
- Default scoring allows gaps

### Local Mode
Soft-clips ends to achieve better alignment:
- Best for: Reads with adapters, ChIP-seq, ATAC-seq
- More sensitive to partial matches

```bash
# Compare modes
bowtie2 --end-to-end -x index -U reads.fq -S end_to_end.sam
bowtie2 --local -x index -U reads.fq -S local.sam
```

## Pre-built Indices

Download pre-built indices to save time:
- Illumina iGenomes: https://support.illumina.com/sequencing/sequencing_software/igenome.html
- Bowtie2 index collections: various sources

## Application-Specific Settings

### ChIP-seq
```bash
bowtie2 -p 8 --very-sensitive --no-mixed --no-discordant \
    -x index -1 chip_R1.fq.gz -2 chip_R2.fq.gz 2> stats.txt | \
    samtools view -bS -q 30 -F 1804 - | \
    samtools sort -o chip.bam -
```

Flags -F 1804 removes:
- Unmapped (4)
- Mate unmapped (8)
- Not primary (256)
- Failed QC (512)
- Duplicate (1024)

### ATAC-seq
```bash
bowtie2 -p 8 --very-sensitive -X 2000 --no-mixed --no-discordant \
    -x index -1 atac_R1.fq.gz -2 atac_R2.fq.gz 2> stats.txt | \
    samtools view -bS -q 30 - | \
    samtools sort -o atac.bam -
```

### Metagenomics
```bash
bowtie2 -p 8 --very-sensitive-local -k 10 \
    -x index -1 r1.fq.gz -2 r2.fq.gz -S aligned.sam 2> stats.txt
```

## Handling Multi-mappers

```bash
# Default: report best alignment only
bowtie2 -x index -U reads.fq -S aligned.sam

# Report up to 5 alignments per read
bowtie2 -k 5 -x index -U reads.fq -S aligned.sam

# Report all alignments (can be large)
bowtie2 -a -x index -U reads.fq -S aligned.sam
```

## Memory and Performance

- Index size: ~3-4GB for human genome
- Memory usage: ~4GB minimum
- Speed scales linearly with threads up to ~16

```bash
# Optimize for speed
bowtie2 --very-fast -p 16 -x index -1 r1.fq -2 r2.fq -S aligned.sam

# Optimize for sensitivity
bowtie2 --very-sensitive -p 8 -x index -1 r1.fq -2 r2.fq -S aligned.sam
```

## Troubleshooting

### Low Alignment Rate
1. Check read quality with FastQC
2. Try local mode: `--local`
3. Increase sensitivity: `--very-sensitive`
4. Verify correct reference genome

### Wrong Fragment Sizes
```bash
# Adjust -X for your library
bowtie2 -X 1000 -x index -1 r1.fq -2 r2.fq -S aligned.sam
```

### Index Errors
```bash
# Rebuild index
rm bt2_index.*
bowtie2-build --threads 8 reference.fa bt2_index
```
