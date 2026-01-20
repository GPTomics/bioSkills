# Bowtie2 Alignment - Usage Guide

## Overview
Bowtie2 is a fast and memory-efficient aligner for short reads. It supports both end-to-end and local alignment modes, making it versatile for ChIP-seq, ATAC-seq, and general short-read alignment.

## Prerequisites
```bash
conda install -c bioconda bowtie2 samtools
```

## Quick Start
Tell your AI agent what you want to do:
- "Align my ChIP-seq reads using Bowtie2"
- "Build a Bowtie2 index for my reference genome"
- "Align ATAC-seq data with appropriate settings"

## Example Prompts

### Basic Alignment
> "Align r1.fq.gz and r2.fq.gz to my reference using Bowtie2"

> "Build a Bowtie2 index for reference.fa"

### ChIP-seq
> "Align my ChIP-seq paired-end reads with high sensitivity and filter low-quality alignments"

> "Run Bowtie2 alignment for ChIP-seq with MAPQ filtering"

### ATAC-seq
> "Align ATAC-seq reads allowing for large fragment sizes up to 2kb"

> "Run Bowtie2 on my ATAC-seq data with appropriate settings for open chromatin"

### Multi-mapping
> "Align reads and report up to 5 alignments per read for repetitive regions"

> "Run Bowtie2 in local mode for reads that may have adapter contamination"

## What the Agent Will Do
1. Build index if not present
2. Select appropriate alignment mode (end-to-end vs local)
3. Apply application-specific parameters (ChIP-seq, ATAC-seq, etc.)
4. Sort and index output BAM file

## Tips
- Use local mode (--local) for reads with potential adapter contamination
- For ChIP-seq/ATAC-seq, use --very-sensitive and filter with MAPQ >= 30
- Set -X appropriately for your expected fragment size distribution
- Use --no-mixed --no-discordant to keep only proper pairs

## When to Use Bowtie2 vs BWA

| Use Case | Bowtie2 | BWA-MEM2 |
|----------|---------|----------|
| ChIP-seq | Preferred | Works |
| ATAC-seq | Preferred | Works |
| RNA-seq | No (use STAR) | No |
| WGS | Works | Preferred |
| Local alignment | Yes | Limited |

## Application-Specific Settings

### ChIP-seq
```bash
bowtie2 -p 8 --very-sensitive --no-mixed --no-discordant \
    -x index -1 chip_R1.fq.gz -2 chip_R2.fq.gz 2> stats.txt | \
    samtools view -bS -q 30 -F 1804 - | \
    samtools sort -o chip.bam -
```

### ATAC-seq
```bash
bowtie2 -p 8 --very-sensitive -X 2000 --no-mixed --no-discordant \
    -x index -1 atac_R1.fq.gz -2 atac_R2.fq.gz 2> stats.txt | \
    samtools view -bS -q 30 - | \
    samtools sort -o atac.bam -
```

## Troubleshooting

### Low Alignment Rate
- Check read quality with FastQC
- Try local mode: --local
- Increase sensitivity: --very-sensitive
- Verify correct reference genome

### Wrong Fragment Sizes
```bash
bowtie2 -X 1000 -x index -1 r1.fq -2 r2.fq -S aligned.sam
```
