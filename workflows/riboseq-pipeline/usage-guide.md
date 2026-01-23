# Ribo-seq Pipeline - Usage Guide

## Overview

Complete workflow from Ribo-seq FASTQ to translation efficiency and ORF detection.

## Prerequisites

```bash
# CLI tools
conda install -c bioconda cutadapt bowtie2 star sortmerna

# Python
pip install plastid ribocode

# R
BiocManager::install('riborex')
```

## Quick Start

- "Analyze my Ribo-seq data from FASTQ to translation efficiency"
- "Run the ribosome profiling pipeline"
- "Detect translated ORFs from my Ribo-seq"

## Example Prompts

### Full Pipeline

> "Run the complete Ribo-seq pipeline"

> "Calculate translation efficiency from my ribosome profiling data"

### Specific Steps

> "Just run P-site calibration"

> "Detect novel ORFs from my Ribo-seq data"

## What the Agent Will Do

1. Trim adapters and filter by size
2. Remove rRNA contamination
3. Align to transcriptome
4. Calibrate P-site offsets
5. Calculate translation efficiency
6. Call translated ORFs
7. Generate QC reports

## Tips

- **Read length** - 28-30nt typical for ribosome-protected fragments
- **3-nt periodicity** - Key QC metric; should see strong triplet periodicity
- **P-site offset** - Usually 12-13nt from 5' end for 28-30nt reads
- **rRNA removal** - Critical; >80% can be rRNA in raw data
- **Paired RNA-seq** - Required for translation efficiency calculation
