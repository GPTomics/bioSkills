# Peak Calling - Usage Guide

## Overview

Call peaks from ChIP-seq data using MACS3 to identify transcription factor binding sites or histone modification regions.

## Prerequisites

```bash
# Conda (recommended)
conda install -c bioconda macs3

# Pip
pip install macs3
```

## Quick Start

Tell your AI agent what you want to do:
- "Call peaks from my ChIP-seq BAM file"
- "Identify H3K27ac enriched regions"
- "Find transcription factor binding sites with input control"

## Example Prompts

### Narrow Peaks (TFs, H3K4me3)
> "Call narrow peaks for my transcription factor ChIP-seq"

> "Run MACS3 on my H3K27ac ChIP-seq with input control"

### Broad Peaks (H3K27me3, H3K36me3)
> "Call broad peaks for my H3K27me3 ChIP-seq data"

> "Identify H3K9me3 domains using broad peak calling"

### ATAC-seq
> "Call peaks from my ATAC-seq data"

### Troubleshooting
> "My peak calling found no peaks, help me troubleshoot"

> "I'm getting too many peaks, how can I increase stringency?"

## What the Agent Will Do

1. Sort and index BAM files if not already done
2. Select appropriate peak type (narrow vs broad) based on target
3. Run MACS3 with appropriate parameters and genome size
4. Generate narrowPeak/broadPeak files and summit BED files
5. Report peak counts and quality metrics
6. Suggest parameter adjustments if results are unexpected

## Tips

- Always use input/IgG control for TF ChIP-seq
- Use `--broad` for histone marks that form domains (H3K27me3, H3K36me3)
- Default q-value threshold is 0.05; use `-q 0.01` for higher stringency
- Use `--nomodel --extsize 200` if model building fails on small datasets
- Check peak numbers: TFs typically have 1,000-50,000 peaks; broad marks can have fewer but larger regions
