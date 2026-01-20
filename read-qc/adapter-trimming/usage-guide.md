# Adapter Trimming - Usage Guide

## Overview
Sequencing adapters must be removed before alignment to prevent misalignment and artifacts. Cutadapt offers precise control with error-tolerant matching, while Trimmomatic provides efficient paired-end handling with palindrome mode.

## Prerequisites
```bash
# Cutadapt
pip install cutadapt
# or
conda install -c bioconda cutadapt

# Trimmomatic
conda install -c bioconda trimmomatic
```

## Quick Start
Tell your AI agent what you want to do:
- "Remove adapters from my paired-end reads"
- "Trim Illumina TruSeq adapters from my FASTQ files"
- "My FastQC shows adapter contamination, please fix it"

## Example Prompts

### Standard Trimming
> "Remove TruSeq adapters from my paired-end FASTQ files using Cutadapt"

> "Trim Nextera adapters from my reads and keep only reads longer than 20bp"

### Troubleshooting
> "FastQC still shows adapters after trimming, increase the error tolerance"

> "I don't know which adapters were used, identify and remove them"

### Tool Selection
> "Use Trimmomatic with palindrome mode for better adapter detection on my overlapping reads"

> "Compare Cutadapt vs Trimmomatic output for my dataset"

## What the Agent Will Do
1. Identify the adapter sequences (from FastQC or library prep info)
2. Select appropriate tool (Cutadapt for most cases, Trimmomatic for palindrome mode)
3. Run trimming with proper parameters for paired-end data
4. Set minimum read length to discard too-short fragments
5. Verify adapter removal with post-trimming QC

## Common Adapter Sequences

| Library | Adapter |
|---------|---------|
| Illumina Universal | `AGATCGGAAGAGC` |
| Nextera | `CTGTCTCTTATACACATCT` |
| Small RNA | `TGGAATTCTCGGGTGCCAAGG` |

## Tips
- Check FastQC "Overrepresented sequences" to identify unknown adapters
- Use `-m 20` or higher to discard very short trimmed reads
- Trimmomatic palindrome mode works best for paired reads with adapter read-through
- Multiple adapters can be specified in a single Cutadapt command
- Always run FastQC after trimming to verify adapter removal

## Resources
- [Cutadapt Documentation](https://cutadapt.readthedocs.io/)
- [Trimmomatic Manual](http://www.usadellab.org/cms/?page=trimmomatic)
- [Illumina Adapter Sequences](https://support.illumina.com/downloads/illumina-adapter-sequences-document-1000000002694.html)
