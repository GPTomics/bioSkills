# Proximity Operations - Usage Guide

## Overview
Proximity operations help you find relationships between genomic features based on their distance. This skill covers finding nearest features (closest), searching within distance windows (window), and extending intervals (slop, flank).

## Prerequisites
```bash
# bedtools (required)
conda install -c bioconda bedtools

# pybedtools (optional, for Python)
pip install pybedtools

# samtools (for creating genome file from FASTA)
conda install -c bioconda samtools
```

## Quick Start
Tell your AI agent what you want to do:
- "Find the nearest gene to each of my ChIP-seq peaks"
- "Get all peaks within 10kb of a TSS"
- "Create promoter regions 2kb upstream of each gene"

## Example Prompts

### Finding Nearest Features
> "Find the nearest gene to each peak in my ChIP-seq data"
> "Assign each enhancer to its closest TSS"
> "Report the distance from each variant to the nearest exon"

### Distance-Based Filtering
> "Get all peaks within 5kb of any promoter"
> "Find variants that are within 1kb of a splice site"
> "Filter enhancers to keep only those within 100kb of a gene"

### Extending Intervals
> "Extend my peak summits by 250bp on each side"
> "Create 2kb promoter regions upstream of each TSS"
> "Get 500bp flanking regions on both sides of my peaks"

### Enhancer-Gene Pairing
> "Pair each enhancer with genes within 100kb"
> "Create enhancer-promoter pairs based on proximity"

## What the Agent Will Do
1. Prepare genome file if needed (chromosome sizes)
2. Extract or create reference features (e.g., TSS from genes)
3. Run appropriate proximity operation (closest, window, slop)
4. Filter results by distance if requested
5. Output paired features with distances

## Choosing the Right Operation

| Goal | Operation | Example |
|------|-----------|---------|
| Find nearest feature | closest | Assign peaks to genes |
| Features within distance | window | Peaks within 10kb of TSS |
| Extend interval borders | slop | Expand peaks by 100bp |
| Get flanking regions | flank | Get regions beside features |
| Move intervals | shift | Shift by fixed distance |

## Tips
- Most operations require a genome file (chromosome sizes); create with `cut -f1,2 reference.fa.fai > genome.txt`
- Use `-s` flag with slop/flank to respect strand direction
- The closest operation returns "." when no feature exists on the same chromosome
- Use `-d` flag with closest to report distances
- For window operations, use `-w` for symmetric windows or `-l`/`-r` for asymmetric

## Resources
- [bedtools closest](https://bedtools.readthedocs.io/en/latest/content/tools/closest.html)
- [bedtools window](https://bedtools.readthedocs.io/en/latest/content/tools/window.html)
- [bedtools slop](https://bedtools.readthedocs.io/en/latest/content/tools/slop.html)
- [bedtools flank](https://bedtools.readthedocs.io/en/latest/content/tools/flank.html)
