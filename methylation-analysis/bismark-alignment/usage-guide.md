# Bismark Alignment - Usage Guide

## Overview
Bismark is the standard tool for aligning bisulfite sequencing reads, handling the complexity of bisulfite conversion (C->T) by creating in-silico converted genomes and performing four-way alignment.

## Prerequisites
```bash
conda install -c bioconda bismark bowtie2 samtools trim-galore

bismark --version
bowtie2 --version
```

## Quick Start
Tell your AI agent what you want to do:
- "Align my WGBS FASTQ files to the hg38 reference genome"
- "Prepare a bisulfite genome index for my reference"
- "Run Bismark alignment with deduplication for my paired-end samples"

## Example Prompts
### Genome Preparation
> "Create a Bismark genome index from my hg38.fa reference file"

> "Prepare a bisulfite-converted genome for mm10"

### Read Alignment
> "Align my paired-end bisulfite sequencing reads to the genome"

> "Run Bismark on my RRBS data without deduplication"

> "Align my non-directional bisulfite library to hg38"

### Troubleshooting
> "My Bismark mapping rate is very low, help me diagnose the issue"

> "Check if my library is directional or non-directional"

## What the Agent Will Do
1. Prepare the bisulfite genome index (if not already done)
2. Trim adapter sequences with Trim Galore
3. Align reads with Bismark using appropriate parameters for your library type
4. Deduplicate aligned reads (for WGBS, skip for RRBS)
5. Generate alignment reports and QC statistics

## Output Files

| File | Description |
|------|-------------|
| *_bismark_bt2.bam | Aligned reads with XM tag (methylation call) |
| *_SE_report.txt | Alignment statistics |
| *.deduplicated.bam | After deduplication |

## Tips
- Always run genome preparation once before alignment (bismark_genome_preparation --bowtie2 genome/)
- Use Trim Galore for adapter trimming before alignment
- Skip deduplication for RRBS data (PCR duplicates expected at restriction sites)
- Check bisulfite conversion efficiency in the alignment report (should be >99%)
- For non-directional libraries, add --non_directional flag
- For PBAT libraries, use --pbat flag
- If mapping rate is low, try relaxing alignment with -N 1
- Use --parallel for multi-core processing instead of running multiple instances
