# m6Anet Analysis - Usage Guide

## Overview

Detect m6A modifications directly from Oxford Nanopore direct RNA sequencing signal.

## Prerequisites

```bash
pip install m6anet ont-fast5-api
conda install -c bioconda minimap2 samtools
```

## Quick Start

- "Detect m6A from my Nanopore direct RNA data"
- "Run m6Anet on my FAST5 files"
- "Find m6A modifications without MeRIP"

## Example Prompts

### Basic Analysis

> "Preprocess my Nanopore FAST5 files for m6Anet"

> "Run m6A detection on direct RNA-seq data"

### Filtering Results

> "Filter m6Anet results for high-confidence sites"

> "Compare m6A between samples from Nanopore data"

## What the Agent Will Do

1. Verify FAST5 files and reference transcriptome
2. Run m6Anet preprocessing (feature extraction)
3. Perform m6A probability inference
4. Filter high-confidence sites
5. Export results with genomic coordinates

## Tips

- **Direct RNA** - Must use ONT direct RNA protocol (not cDNA)
- **FAST5 required** - Needs raw signal data, not just FASTQ
- **Transcriptome ref** - Use transcript sequences, not genome
- **Coverage** - Higher coverage improves detection accuracy
- **Probability > 0.9** - Conservative threshold for high-confidence calls
