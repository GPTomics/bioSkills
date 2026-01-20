# Structural Variant Detection - Usage Guide

## Overview
Long reads excel at detecting structural variants (deletions, insertions, inversions, duplications) that are difficult or impossible to detect with short reads. Sniffles2 is the most popular tool, with cuteSV and SVIM as alternatives.

## Prerequisites
```bash
conda install -c bioconda sniffles cutesv svim pbsv bcftools
```

## Quick Start
Tell your AI agent what you want to do:
- "Call structural variants from my ONT alignment"
- "Detect large insertions and deletions using Sniffles"

## Example Prompts

### Basic SV Calling
> "Call SVs from aligned.bam using Sniffles2 and output to structural_variants.vcf"

> "Run cuteSV on my PacBio HiFi alignment for SV detection"

### With Reference for Insertions
> "Call SVs with Sniffles and include the reference so insertion sequences are output"

### Filtering
> "Call SVs and filter for high-confidence variants with QUAL >= 20"

### Population Calling
> "Merge SV calls from multiple samples for population analysis"

### Tandem Repeat Regions
> "Call SVs with tandem repeat annotation for improved accuracy in repetitive regions"

## What the Agent Will Do
1. Verify alignment file is sorted and indexed
2. Run SV caller with appropriate settings for your platform (ONT vs HiFi)
3. Apply quality filtering
4. Optionally merge calls from multiple samples

## Tool Comparison

| Tool | Best For | Notes |
|------|----------|-------|
| Sniffles2 | General purpose | Most popular, population calling |
| cuteSV | High sensitivity | Platform-specific settings |
| SVIM | ONT data | Good insertion detection |
| pbsv | PacBio HiFi | PacBio-optimized |

## Tips
- Need at least 10x coverage for reliable SV calling
- Always provide reference with `--reference` to get insertion sequences in output
- Use tandem repeat BED file with `--tandem-repeats` for better accuracy in repetitive regions
- Filter by QUAL >= 20 for high-confidence calls
- For PacBio HiFi, consider pbsv which is optimized for that platform

## Resources
- [Sniffles2 GitHub](https://github.com/fritzsedlazeck/Sniffles)
- [cuteSV GitHub](https://github.com/tjiangHIT/cuteSV)
