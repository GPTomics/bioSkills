# TCR/BCR Pipeline - Usage Guide

## Overview

Complete workflow from TCR/BCR sequencing FASTQ files to clonotype analysis and diversity metrics.

## Prerequisites

```bash
# MiXCR (requires license)
# Download from https://github.com/milaboratory/mixcr

# VDJtools
conda install -c bioconda vdjtools
```

## Quick Start

- "Analyze my TCR-seq data from FASTQ to diversity"
- "Run the complete immune repertoire pipeline"
- "Process bulk TCR sequencing end-to-end"

## Example Prompts

### Full Pipeline

> "Run the TCR pipeline on my paired-end FASTQ files"

> "Analyze my BCR repertoire from sequencing to clonotypes"

### Specific Steps

> "Just run MiXCR alignment and assembly"

> "Calculate diversity metrics from my clonotype table"

## What the Agent Will Do

1. Align reads to V(D)J reference (MiXCR)
2. Assemble clonotypes
3. Export clonotype table
4. Calculate diversity metrics (VDJtools)
5. Generate visualization plots
6. Report QC metrics at each step

## Tips

- **Species** - Use -s hsa (human) or -s mmu (mouse)
- **Protocol** - rna-seq for bulk, kTCR-seq for targeted
- **Clonotype definition** - Default is CDR3+V+J
- **Diversity** - Shannon entropy, Chao1, clonality metrics
- **Overlap** - F2 metric for public clone detection
