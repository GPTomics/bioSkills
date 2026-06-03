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

1. **Alignment** - Map Hi-C read pairs with bwa-mem2 -SP5M (mates aligned independently)
2. **Pairs** - Classify, filter, and deduplicate; judge library quality from long-range cis
3. **Matrix** - Build a cooler and zoomify to a multi-resolution .mcool
4. **Balance** - ICE normalization (required before any analysis)
5. **Compartments** - A/B eigenvector at 100kb, sign-phased by GC
6. **TADs** - Insulation-score boundaries across a window sweep at 10kb
7. **Loops** - Dot calling at 10kb IF the map is deep enough; else APA on known anchors

## Tips

- **Depth dictates the feature**: compartments are cheap, TAD boundaries need a moderate map, de-novo loops need billions of contacts (~5B in Rao 2014). Decide the resolution from the depth, not the other way around.
- **Resolution to feature**: 100kb-1Mb for compartments, 10-40kb for TADs, 5-10kb for loops (1-2kb for Micro-C).
- **Phasing is not optional**: the compartment eigenvector sign is arbitrary until oriented by a GC or gene-density track.
- **QC**: the long-range cis (>=20kb) fraction is the one-number library readout; trans is a noise floor whose threshold is genome-size-dependent. High duplicate rate = low complexity (not fixable by sequencing deeper).
- **Protein-directed assays branch out**: HiChIP, PLAC-seq, and Capture Hi-C need FitHiChIP/MAPS/CHiCAGO (hi-c-analysis/hichip-plac-loops), not generic dot calling.
