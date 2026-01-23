# FASTQ Quality - Usage Guide

## Overview

This skill enables AI agents to help you work with FASTQ quality scores for NGS data analysis using Biopython.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Calculate average quality for each read"
- "Filter reads with mean quality below 25"
- "Trim low-quality bases from read ends"
- "Show quality statistics for my FASTQ file"

## Example Prompts

### Quality Analysis
> "Show me the quality distribution of reads.fastq"

### Filtering
> "Keep only reads with average quality >= 30"

### Trimming
> "Trim bases from the 3' end where quality drops below 20"

### Statistics
> "Generate per-position quality profile for the first 50 bases"

## What the Agent Will Do
1. Parse FASTQ records including quality scores
2. Decode quality encoding (Phred33/64)
3. Calculate per-base or per-read quality metrics
4. Report summary statistics

## Quality Score Reference

- Q20 = 99% accuracy (1 error per 100 bases)
- Q30 = 99.9% accuracy (1 error per 1000 bases)
- Q40 = 99.99% accuracy (1 error per 10000 bases)

## Tips

- Most modern FASTQ uses Phred+33 encoding (format: 'fastq')
- For old Illumina data, try 'fastq-illumina' format
- Process large files as iterators to save memory
- Quality often drops toward 3' end of reads
