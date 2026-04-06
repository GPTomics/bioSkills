# Long-Read SV Pipeline - Usage Guide

## Overview

This workflow detects structural variants (deletions, insertions, inversions, duplications) from Oxford Nanopore or PacBio long-read sequencing data.

## Prerequisites

```bash
conda install -c bioconda minimap2 samtools sniffles cutesv nanoplot bcftools
```

## Quick Start

Tell your AI agent what you want to do:
- "Detect structural variants from my Nanopore data"
- "Run the long-read SV pipeline on my PacBio HiFi reads"
- "Find deletions and insertions in my ONT sequencing"

## Example Prompts

### SV calling
> "Call SVs from my aligned long reads"

> "Use Sniffles to detect structural variants"

> "Find large deletions in my sample"

### Multi-sample
> "Merge SVs across multiple samples"

> "Joint call SVs from my cohort"

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| Long reads | FASTQ | ONT or PacBio reads |
| Reference | FASTA | Reference genome |
| Coverage | >10x | Higher is better for SVs |

## What the Agent Will Do

1. Assess read length and quality with NanoPlot
2. Align reads to reference with minimap2 (map-ont or map-hifi preset)
3. Call structural variants with Sniffles2 or cuteSV
4. Filter SV calls by quality and size
5. Annotate with AnnotSV for gene and clinical context

## Sniffles vs cuteSV

| Feature | Sniffles2 | cuteSV |
|---------|-----------|--------|
| Speed | Moderate | Fast |
| Accuracy | High | High |
| Multi-sample | Built-in | External merge |
| Best for | General use | Large cohorts |

## Tips

- 15-30x coverage recommended for reliable SV calling
- Longer reads detect larger SVs and resolve complex rearrangements better
- Provide tandem repeat annotation to Sniffles2 to improve accuracy
- Start with QUAL>=20, adjust based on validation
- Use PBMM2 instead of minimap2 for PacBio HiFi data with pbsv

## Related Skills

- variant-calling/structural-variant-calling - Short-read SV calling and consensus approach
- long-read-sequencing/structural-variants - Long-read SV details
- long-read-sequencing/long-read-alignment - Minimap2 alignment options
- variant-calling/vcf-manipulation - Merge SV calls across callers
