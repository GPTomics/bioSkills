# long-read-sequencing

## Overview

Analysis of long-read sequencing data from Oxford Nanopore and PacBio. Covers alignment with minimap2, polishing and variant calling with medaka, and structural variant detection with Sniffles.

**Tool type:** cli | **Primary tools:** minimap2, medaka, Sniffles, NanoPlot

## Skills

| Skill | Description |
|-------|-------------|
| long-read-alignment | Align long reads with minimap2 |
| medaka-polishing | Polish assemblies and call variants with medaka |
| structural-variants | Detect SVs from long reads |
| long-read-qc | Quality control for long reads |

## Example Prompts

- "Align my Nanopore reads with minimap2"
- "Map PacBio HiFi reads to the reference genome"
- "Create a sorted BAM from long-read alignment"
- "Polish my assembly with medaka"
- "Call variants from Nanopore reads with medaka"
- "Generate consensus sequence for my region"
- "Find structural variants from my long reads"
- "Detect deletions and insertions with Sniffles"
- "Call SVs from PacBio alignments"
- "Check the quality of my Nanopore reads"
- "Generate read length distribution"
- "Filter reads by quality score"

## Requirements

```bash
# minimap2
conda install -c bioconda minimap2

# medaka
conda install -c bioconda medaka

# SV callers
conda install -c bioconda sniffles cutesv

# QC tools
conda install -c bioconda nanoplot chopper
```

## Related Skills

- **alignment-files** - BAM manipulation
- **variant-calling** - Short-read variant calling
- **sequence-io** - FASTQ handling
