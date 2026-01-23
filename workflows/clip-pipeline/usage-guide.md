# CLIP-seq Pipeline - Usage Guide

## Overview

Complete workflow from CLIP-seq FASTQ to binding sites, annotation, and motif enrichment.

## Prerequisites

```bash
conda install -c bioconda umi_tools cutadapt star clipper bedtools homer
```

## Quick Start

- "Analyze my CLIP-seq data for binding sites"
- "Run the protein-RNA interaction pipeline"
- "Process my eCLIP data end-to-end"

## Example Prompts

### Full Pipeline

> "Run the complete CLIP-seq pipeline"

> "Find RBP binding sites and enriched motifs"

### Specific Steps

> "Just call peaks from my aligned CLIP BAM"

> "Run motif analysis on my binding sites"

## What the Agent Will Do

1. Extract UMIs if present
2. Trim adapters and quality filter
3. Align to genome (STAR)
4. Deduplicate using UMIs
5. Call binding site peaks (CLIPper)
6. Annotate to transcriptomic features
7. Find enriched motifs (HOMER)

## Tips

- **UMI** - Critical for PCR duplicate removal in CLIP
- **Input control** - SMInput recommended for eCLIP
- **Peak width** - CLIP peaks typically narrow (20-50nt)
- **Crosslink sites** - Expect enrichment at specific nucleotides
- **DRACH/RBP motifs** - Known motif should be enriched
- **Replicate overlap** - Use peaks present in 2+ replicates
