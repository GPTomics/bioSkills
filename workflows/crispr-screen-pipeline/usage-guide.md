# CRISPR Screen Pipeline Usage Guide

## Overview

This workflow takes raw FASTQ files from pooled CRISPR screens through guide counting, quality control, statistical analysis, and hit calling to identify genes affecting a phenotype.

## Prerequisites

```bash
pip install mageck mageck-vispr
conda install -c bioconda mageck
```

## Quick Start

Tell your AI agent what you want to do:
- "Run the CRISPR screen pipeline on my FASTQ files"
- "Count guides and run MAGeCK on my dropout screen"
- "Analyze my CRISPRi screen for essential genes"

## Example Prompts

### Basic Analysis
> "I have FASTQ files from a CRISPR knockout screen, run the full pipeline"

> "Process my pooled screen with MAGeCK RRA"

### Screen Types
> "Analyze my drug resistance screen for enriched guides"

> "Find essential genes in my dropout screen using BAGEL2"

### Custom Analysis
> "Run MAGeCK MLE with my multi-condition design"

> "Compare my early vs late timepoints in a fitness screen"

## When to Use This Pipeline

- Pooled CRISPR knockout/knockdown screens
- Dropout (negative selection) screens
- Enrichment (positive selection) screens
- CRISPRi/CRISPRa screens
- Drug resistance screens
- Fitness/essentiality screens

## Required Inputs

1. **FASTQ files** - Sequencing reads containing guide sequences
2. **Library file** - CSV with sgRNA sequences and gene mappings
3. **Sample annotation** - Which samples are controls/treatments

## Library File Format

```csv
sgRNA,Gene,Sequence
sgRNA1,TP53,ACGTACGTACGTACGTACGT
sgRNA2,TP53,GCTAGCTAGCTAGCTAGCTA
sgRNA3,BRCA1,ATCGATCGATCGATCGATCG
```

## Pipeline Steps

### 1. Guide Counting
- Extracts and counts guide sequences from FASTQ
- Outputs count matrix (guides × samples)

### 2. Quality Control
- Check mapping rates (>70%)
- Calculate zero-count guides (<20%)
- Assess distribution evenness (Gini index <0.4)
- Verify replicate correlation

### 3. Statistical Analysis
- **MAGeCK RRA**: Robust rank aggregation, good for most screens
- **MAGeCK MLE**: Maximum likelihood, for complex designs
- **BAGEL2**: Bayesian analysis, good for essentiality screens

### 4. Hit Calling
- FDR threshold (typically <0.05)
- Effect size filter (|LFC| > 0.5)
- Combine multiple methods for confidence

## Screen Types

### Negative Selection (Dropout)
Guides targeting essential genes drop out over time.
- Use `neg|fdr` and `neg|lfc` columns
- Hits have negative LFC

### Positive Selection (Enrichment)
Guides targeting resistance genes become enriched.
- Use `pos|fdr` and `pos|lfc` columns
- Hits have positive LFC

## Common Issues

### Low mapping rate
- Check if reads are trimmed correctly
- Verify library sequences match

### High Gini index
- PCR amplification bias
- Consider resequencing

### No significant hits
- Screen may not have worked
- Adjust thresholds cautiously
- Check positive controls

## Output Files

| File | Description |
|------|-------------|
| experiment.count.txt | Guide count matrix |
| negative_screen.gene_summary.txt | Gene-level statistics |
| negative_screen.sgrna_summary.txt | Guide-level statistics |
| negative_hits.csv | Called hit genes |
| volcano_plot.png | Visualization |

## Tips

- **Library file**: Ensure sgRNA sequences match your library exactly
- **Replicates**: 3+ biological replicates improve hit confidence
- **Normalization**: Use non-targeting controls if available
- **Multiple methods**: Run both MAGeCK and BAGEL2 for confident hits
- **Positive controls**: Include known essential/drug-target genes to validate screen

## References

- MAGeCK: doi:10.1186/s13059-014-0554-4
- MAGeCKFlute: doi:10.1038/s41596-018-0113-7
- BAGEL2: doi:10.1186/s13073-019-0698-6
