# Paired-End FASTQ - Usage Guide

## Overview

This skill enables AI agents to help you work with paired-end sequencing data (R1/R2 files) using Biopython.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Filter paired reads keeping only pairs where both pass quality threshold"
- "Interleave R1 and R2 files into a single file"
- "Deinterleave this combined FASTQ into separate R1/R2 files"
- "Find all paired files in this directory"

## Example Prompts

### Filtering
> "Filter my paired FASTQ files, keeping pairs where both reads have mean quality >= 30"

### Interleaving
> "Combine reads_R1.fastq and reads_R2.fastq into an interleaved file"

### Deinterleaving
> "Split this interleaved FASTQ into separate R1 and R2 files"

### Validation
> "Check if my R1 and R2 files have matching read counts and IDs"

## Common Naming Patterns

| Pattern | R1 | R2 |
|---------|-----|-----|
| Illumina | `sample_R1_001.fastq` | `sample_R2_001.fastq` |
| Simple | `sample_1.fastq` | `sample_2.fastq` |
| Underscore | `sample_R1.fastq` | `sample_R2.fastq` |

## What the Agent Will Do
1. Open both R1 and R2 FASTQ files
2. Iterate through paired records simultaneously
3. Verify read pairing by ID matching
4. Process paired reads together

## Tips

- Always filter pairs together (both must pass)
- Verify pair counts match before processing
- Use streaming for large files to avoid memory issues
- Interleaved format: R1, R2, R1, R2, ... alternating
- Most tools accept either separate or interleaved input
