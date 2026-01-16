# Bismark Alignment - Usage Guide

## Overview

Bismark is the standard tool for aligning bisulfite sequencing reads. It handles the complexity of bisulfite conversion (C→T) by creating in-silico converted genomes and performing four-way alignment.

## When to Use This Skill

- You have raw bisulfite sequencing FASTQ files
- You need to align WGBS, RRBS, or targeted bisulfite-seq data
- You want to prepare data for methylation calling

## Installation

```bash
# Conda (recommended)
conda install -c bioconda bismark bowtie2 samtools

# Verify
bismark --version
bowtie2 --version
```

## Workflow

### 1. Prepare Genome (Once per Reference)

```bash
# Download reference genome
wget https://example.com/hg38.fa.gz
gunzip hg38.fa.gz
mkdir genome && mv hg38.fa genome/

# Create bisulfite genome index
bismark_genome_preparation --bowtie2 genome/
```

### 2. Trim Reads (Recommended)

```bash
# Trim Galore handles adapter trimming
trim_galore --illumina --paired \
    sample_R1.fastq.gz sample_R2.fastq.gz
```

### 3. Align Reads

```bash
bismark --genome genome/ \
    -1 sample_R1_val_1.fq.gz \
    -2 sample_R2_val_2.fq.gz \
    -o aligned/
```

### 4. Deduplicate (WGBS Only)

```bash
deduplicate_bismark --paired --bam aligned/sample_bismark_bt2_pe.bam
```

## Common Issues

### Low Mapping Rate

- Check bisulfite conversion efficiency in report
- Verify correct library type (directional vs non-directional)
- Try relaxing alignment parameters (-N 1)

### Memory Issues

- Use --parallel instead of running multiple instances
- Increase --temp_dir to a disk with more space

### Wrong Library Type

Symptoms: Very low mapping, biased methylation calls
Solution: Try --non_directional or --pbat

## Output Files

| File | Description |
|------|-------------|
| *_bismark_bt2.bam | Aligned reads with XM tag (methylation call) |
| *_SE_report.txt | Alignment statistics |
| *.deduplicated.bam | After deduplication |

## Resources

- [Bismark User Guide](https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html)
- [Bismark GitHub](https://github.com/FelixKrueger/Bismark)
