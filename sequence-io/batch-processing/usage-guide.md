# Batch Processing - Usage Guide

## Overview

This skill enables AI agents to help you process multiple sequence files efficiently using Biopython.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Merge all FASTA files in the data folder"
- "Split this large file into 1000-record chunks"
- "Convert all GenBank files to FASTA"
- "Count sequences in each file in the directory"

## Example Prompts

### Merging
> "Combine all .fasta files in samples/ into one file"

### Splitting
> "Split large.fasta into files of 500 sequences each"

### Batch Conversion
> "Convert all .gb files in genbank/ to FASTA format in fasta/"

### Statistics
> "Generate a summary CSV with sequence counts and lengths for all FASTA files"

### Parallel Processing
> "Process all FASTQ files in parallel and count reads"

## What the Agent Will Do
1. Load sequences from the input files
2. Apply the requested processing to each sequence
3. Track progress across all files
4. Output results to specified location

## Tips

- Use pathlib.Path for cross-platform file handling
- Use generators when merging to avoid loading all into memory
- For CPU-intensive tasks, use multiprocessing.Pool
- For I/O-bound tasks, ThreadPoolExecutor may be faster
- Always use rglob() for recursive directory search
