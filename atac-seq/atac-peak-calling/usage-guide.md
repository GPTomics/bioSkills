# ATAC-seq Peak Calling - Usage Guide

## Overview
Call accessible chromatin regions from ATAC-seq data using MACS3 or Genrich, with specialized handling for Tn5 transposase cut sites and nucleosome-free region detection.

## Prerequisites
```bash
pip install macs3
conda install -c bioconda samtools bedtools genrich
```

## Quick Start
Tell your AI agent what you want to do:
- "Call ATAC-seq peaks from my BAM file"
- "Separate nucleosome-free and mono-nucleosomal reads before peak calling"

## Example Prompts
### Basic Peak Calling
> "Call peaks from my ATAC-seq BAM file using MACS3 with proper Tn5 offset correction"

### Fragment-Based Analysis
> "Separate my ATAC-seq reads by fragment size and call peaks on nucleosome-free regions only"

### Consensus Peaks
> "Create a consensus peak set from multiple ATAC-seq samples"

### Quality Filtering
> "Filter my ATAC-seq peaks by q-value and remove blacklist regions"

## What the Agent Will Do
1. Filter BAM for properly paired, high-quality reads and remove mitochondrial reads
2. Optionally separate reads by fragment size (NFR, mono-nucleosomal)
3. Call peaks with MACS3 using ATAC-seq specific parameters (--shift -75 --extsize 150 --keep-dup all)
4. Filter peaks by quality and remove blacklist regions
5. Generate consensus peak set if multiple samples provided

## Key Differences from ChIP-seq

| Aspect | ATAC-seq | ChIP-seq |
|--------|----------|----------|
| Signal source | Tn5 cut sites | Protein binding |
| Control | No input control | Input/IgG required |
| Fragment pattern | Nucleosome periodicity | Smooth enrichment |
| Duplicates | Keep all | Remove PCR duplicates |

## Fragment Size Selection

| Fragment Size | Origin | Analysis Use |
|---------------|--------|--------------|
| <100 bp | Nucleosome-free | TF binding, footprinting |
| 180-247 bp | Mono-nucleosome | Nucleosome positioning |
| 315-473 bp | Di-nucleosome | Chromatin structure |
| 558-615 bp | Tri-nucleosome | Chromatin structure |

## Tips
- Use `--shift -75 --extsize 150` for paired-end ATAC-seq to account for Tn5 offset
- Always use `--keep-dup all` since ATAC-seq has legitimate duplicates
- Filter for NFR reads (<100bp) when looking for TF binding sites
- Remove ENCODE blacklist regions to reduce false positives
- For consensus peaks across samples, require presence in at least 2 samples
