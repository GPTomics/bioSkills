# Hi-C Pipeline - Usage Guide

## Overview

This workflow analyzes Hi-C chromosome conformation capture data to identify compartments, TADs, and chromatin loops.

## Prerequisites

```bash
conda install -c bioconda bwa-mem2 pairtools cooler
pip install cooltools
```

## Quick Start

Tell your AI agent what you want to do:
- "Analyze my Hi-C data for TADs and loops"
- "Find A/B compartments in my Hi-C matrix"
- "Process Hi-C FASTQ to contact matrix"

## Example Prompts

### Processing
> "Generate a contact matrix from my Hi-C pairs"

> "Balance my cooler file with ICE"

### Analysis
> "Call TADs using insulation score"

> "Detect chromatin loops"

> "Find compartment boundaries"

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| Hi-C reads | FASTQ | Paired-end Hi-C library |
| Reference | FASTA | Genome reference |
| Chromosome sizes | TSV | For cooler |

## What the Workflow Does

1. **Alignment** - Map Hi-C read pairs
2. **Pairs** - Filter and deduplicate
3. **Matrix** - Generate contact matrix
4. **Balance** - ICE normalization
5. **Compartments** - A/B eigenvector
6. **TADs** - Insulation score boundaries
7. **Loops** - Dot calling

## Tips

- **Resolution**: 100kb for compartments, 10kb for TADs, 5-10kb for loops
- **Sequencing depth**: 500M-1B reads for comprehensive analysis
- **QC**: Check cis/trans ratio and duplicate rate
