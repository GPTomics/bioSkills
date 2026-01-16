# Peak Calling with MACS3 - Usage Guide

## Overview

MACS3 (Model-based Analysis of ChIP-Seq) is the standard tool for identifying enriched regions (peaks) in ChIP-seq data. It models the fragment size, uses local background correction, and provides statistical significance for each peak.

## When to Use This Skill

- You have aligned ChIP-seq BAM files
- You need to identify binding sites or enriched regions
- You want narrowPeak or broadPeak files for downstream analysis

## Installation

```bash
# Conda (recommended)
conda install -c bioconda macs3

# Pip
pip install macs3

# Verify
macs3 --version
```

## Choosing Peak Type

| Target | Peak Type | MACS3 Flag |
|--------|-----------|------------|
| Transcription factors | Narrow | (default) |
| H3K4me3, H3K27ac | Narrow | (default) |
| H3K4me1 | Narrow or Broad | --broad optional |
| H3K36me3 | Broad | --broad |
| H3K27me3, H3K9me3 | Broad | --broad |
| ATAC-seq | Narrow | --nomodel --extsize 200 |

## Basic Workflow

### 1. Prepare BAM Files

```bash
# Sort and index (if not already done)
samtools sort chip.bam -o chip.sorted.bam
samtools index chip.sorted.bam
```

### 2. Call Peaks

```bash
# Narrow peaks with control
macs3 callpeak -t chip.sorted.bam -c input.sorted.bam \
    -f BAM -g hs -n experiment --outdir peaks/

# Broad peaks
macs3 callpeak -t chip.sorted.bam -c input.sorted.bam \
    -f BAM -g hs -n experiment --broad --outdir peaks/
```

### 3. Check Results

```bash
# Number of peaks
wc -l peaks/experiment_peaks.narrowPeak

# Top peaks by significance
sort -k8,8nr peaks/experiment_peaks.narrowPeak | head
```

## Common Issues

### No Peaks Called

- Check BAM file has reads: `samtools view -c chip.bam`
- Lower q-value threshold: `-q 0.1`
- Check genome size parameter matches your data

### Model Building Failed

- Use `--nomodel --extsize 200` for small datasets
- Use `--nomodel --extsize 147` for MNase-seq

### Too Many Peaks

- Increase q-value stringency: `-q 0.01`
- Use input control (essential for TF ChIP-seq)

## Output Interpretation

### narrowPeak columns

1. chrom
2. chromStart
3. chromEnd
4. name
5. score (0-1000)
6. strand
7. signalValue (fold enrichment)
8. pValue (-log10)
9. qValue (-log10)
10. peak (summit position relative to start)

## Resources

- [MACS3 GitHub](https://github.com/macs3-project/MACS)
- [MACS3 Documentation](https://macs3-project.github.io/MACS/)
