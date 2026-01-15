# SRA Data - Usage Guide

## Overview

This skill enables AI agents to help you download raw sequencing data (FASTQ) from the NCBI Sequence Read Archive using the SRA toolkit command-line tools.

## Prerequisites

Install the SRA toolkit:

```bash
# macOS
brew install sratoolkit

# Ubuntu/Debian
sudo apt install sra-toolkit

# conda (recommended)
conda install -c bioconda sra-tools
```

Verify installation:
```bash
fasterq-dump --version
```

## Quick Start

Tell your AI agent what you want to do:

- "Download the FASTQ files for SRR12345678"
- "Download all runs from BioProject PRJNA123456"
- "Prefetch these SRA accessions and convert to FASTQ"
- "How do I find the SRA accessions for a specific study?"

## Example Prompts

### Single Download
> "Download SRR12345678 as FASTQ files"

### Multiple Downloads
> "Download these SRA runs: SRR1234567, SRR1234568, SRR1234569"

### Large Downloads
> "I need to download 50GB worth of SRA data - what's the best approach?"

### Finding Accessions
> "Find all SRA runs for human RNA-seq from project PRJNA123456"

### Paired-End Data
> "Download paired-end FASTQ files for SRR12345678"

## What the Agent Will Do

1. Use fasterq-dump for standard downloads
2. Recommend prefetch for large files
3. Set appropriate thread counts for your system
4. Handle paired vs single-end data automatically
5. Provide validation commands for integrity checking

## Tips

- Use prefetch first for files >20GB
- fasterq-dump auto-splits paired-end reads
- Default output location is current directory
- SRR* accessions are the download unit (not SRX, SRS, or SRP)
- Configure vdb-config for custom cache locations
