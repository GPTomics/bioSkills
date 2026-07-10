# CNV Pipeline - Usage Guide

## Overview

This workflow detects copy number variants and forks by question: germline CNVs -> GATK gCNV; somatic exome/panel -> CNVkit; WGS -> ASCAT/FACETS/PURPLE; cfDNA -> ichorCNA. It commits the build + target BED + assay-matched panel-of-normals, builds the reference before segmentation, fits purity/ploidy before integer calls, and diploid-centers before GISTIC2. CNVkit (below) is the somatic exome/panel path.

## Prerequisites

```bash
conda install -c bioconda cnvkit

# Or pip
pip install cnvkit
```

## Quick Start

Tell your AI agent what you want to do:
- "Detect CNVs from my exome sequencing data"
- "Run CNVkit on my tumor-normal pair"
- "Find copy number changes in my samples"

## Example Prompts

### CNV calling
> "Call CNVs using my pool of normals"

> "Detect amplifications and deletions"

> "Create a copy number profile for my sample"

### Visualization
> "Plot the CNV scatter for chromosome 17"

> "Create a heatmap of CNVs across samples"

> "Generate a diagram showing CNV locations"

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| BAM files | Aligned reads | Tumor and/or normal |
| Target BED | BED file | Capture regions |
| Reference | FASTA | Genome reference |

## What the Workflow Does

1. **Target Preparation** - Define regions to analyze
2. **Coverage** - Calculate read depth per region
3. **Reference** - Create baseline from normals
4. **Calling** - Identify gains and losses
5. **Visualization** - Plot profiles
6. **Annotation** - Map to genes

## Tips

- **Normals**: Pool of normals improves accuracy
- **Matched normal**: Best for somatic CNV detection
- **Coverage**: Uniform coverage is critical
- **Thresholds**: Adjust based on noise level
