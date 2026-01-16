# Proximity Operations - Usage Guide

## Overview

Proximity operations help you find relationships between genomic features based on their distance. This skill covers finding nearest features (closest), searching within distance windows (window), and extending intervals (slop, flank).

## When to Use This Skill

- Find the nearest gene to each ChIP-seq peak
- Identify peaks within a certain distance of TSS
- Create promoter regions around gene starts
- Extend peak summits to defined window sizes
- Assign regulatory elements to target genes

## Installation

```bash
# bedtools (required)
conda install -c bioconda bedtools

# pybedtools (optional, for Python)
pip install pybedtools
```

## Choosing the Right Operation

| Goal | Operation | Example |
|------|-----------|---------|
| Find nearest feature | closest | Assign peaks to genes |
| Features within distance | window | Peaks within 10kb of TSS |
| Extend interval borders | slop | Expand peaks by 100bp |
| Get flanking regions | flank | Get regions beside features |
| Move intervals | shift | Shift by fixed distance |

## Key Concept: Genome File

Many operations require a genome file specifying chromosome sizes:

```bash
# Create from FASTA index
samtools faidx reference.fa
cut -f1,2 reference.fa.fai > genome.txt

# Or download from UCSC
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes -O genome.txt
```

Format:
```
chr1	248956422
chr2	242193529
...
```

## Basic Workflow

### 1. Prepare TSS File

```bash
# Extract TSS from gene BED file
awk -v OFS='\t' '{
    if ($6 == "+") print $1, $2, $2+1, $4, $5, $6;
    else print $1, $3-1, $3, $4, $5, $6;
}' genes.bed > tss.bed
```

### 2. Find Nearest Gene to Peaks

```bash
bedtools closest -a peaks.bed -b tss.bed -d > peaks_nearest_gene.bed
```

### 3. Filter by Distance

```bash
# Keep only peaks within 10kb of a gene
awk '$NF <= 10000' peaks_nearest_gene.bed > peaks_within_10kb.bed
```

### 4. Create Promoter Regions

```bash
# 2kb upstream of TSS
bedtools slop -i tss.bed -g genome.txt -l 2000 -r 500 -s > promoters.bed
```

## Common Patterns

### Enhancer-Gene Pairing

```python
import pybedtools

enhancers = pybedtools.BedTool('enhancers.bed')
tss = pybedtools.BedTool('tss.bed')

# Find nearest TSS to each enhancer
pairs = enhancers.closest(tss, d=True)

# Filter to within 100kb
nearby = pairs.filter(lambda x: abs(int(x.fields[-1])) <= 100000)
nearby.saveas('enhancer_gene_pairs.bed')
```

### Peaks Near TSS

```bash
bedtools window -a peaks.bed -b tss.bed -w 5000 > peaks_near_tss.bed
```

## Common Issues

### Missing Genome File

**Problem:** "Could not open genome file"
**Solution:** Create genome file from FASTA or download from UCSC

### Coordinates Exceed Chromosome

**Problem:** slop/flank creates negative coordinates
**Solution:** bedtools automatically clips to chromosome boundaries

### No Closest Found

**Problem:** "." returned for closest
**Solution:** No feature on same chromosome; use -D ref for cross-chromosome

## Resources

- [bedtools closest](https://bedtools.readthedocs.io/en/latest/content/tools/closest.html)
- [bedtools window](https://bedtools.readthedocs.io/en/latest/content/tools/window.html)
- [bedtools slop](https://bedtools.readthedocs.io/en/latest/content/tools/slop.html)
- [bedtools flank](https://bedtools.readthedocs.io/en/latest/content/tools/flank.html)
