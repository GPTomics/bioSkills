# BigWig Tracks - Usage Guide

## Overview

BigWig is an indexed binary format for continuous genomic data like coverage, ChIP-seq signal, or conservation scores. It's efficient for genome browsers and programmatic access. This skill covers creating bigWig files from bedGraph and extracting values with pyBigWig.

## When to Use This Skill

- Convert bedGraph to bigWig for browser visualization
- Extract signal values from bigWig for analysis
- Create normalized coverage tracks
- Compare signal between samples

## Installation

```bash
# UCSC tools (CLI conversion)
conda install -c bioconda ucsc-bedgraphtobigwig ucsc-bigwigtobedgraph

# pyBigWig (Python read/write)
pip install pyBigWig

# deepTools (advanced bigWig operations)
conda install -c bioconda deeptools
```

## Key Concepts

### Why BigWig over bedGraph?

| Feature | bedGraph | bigWig |
|---------|----------|--------|
| File size | Large | ~10x smaller |
| Random access | Sequential only | Indexed |
| Browser support | Basic | Full |
| Query speed | Slow | Fast |

### Chromosome Sizes File

Most bigWig operations require a chromosome sizes file:

```bash
# Format: chr<TAB>size
chr1	248956422
chr2	242193529
```

Create from:
```bash
# FASTA index
cut -f1,2 reference.fa.fai > chrom.sizes

# Download from UCSC
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
```

## Basic Workflow

### 1. Generate bedGraph

```bash
bedtools genomecov -ibam alignments.bam -bg > coverage.bedGraph
```

### 2. Sort bedGraph

```bash
sort -k1,1 -k2,2n coverage.bedGraph > coverage.sorted.bedGraph
```

### 3. Convert to bigWig

```bash
bedGraphToBigWig coverage.sorted.bedGraph hg38.chrom.sizes coverage.bw
```

### 4. View in Browser

Upload `coverage.bw` to UCSC Genome Browser or IGV.

## Common Patterns

### Quick Normalized Track

```bash
# Using deepTools (recommended)
bamCoverage -b alignments.bam -o coverage.bw --normalizeUsing CPM
```

### Extract Signal for Regions

```python
import pyBigWig

bw = pyBigWig.open('signal.bw')
mean_signal = bw.stats('chr1', 1000000, 1100000, type='mean')[0]
print(f'Mean signal: {mean_signal}')
bw.close()
```

## Common Issues

### Unsorted bedGraph

**Problem:** "bedGraph not sorted"
**Solution:** Sort by chromosome and position

```bash
sort -k1,1 -k2,2n input.bedGraph > sorted.bedGraph
```

### Chromosome Name Mismatch

**Problem:** Chromosome names in bedGraph don't match chrom.sizes
**Solution:** Ensure consistent naming (chr1 vs 1)

```bash
# Add 'chr' prefix
sed 's/^/chr/' input.bedGraph > with_chr.bedGraph
```

### Missing Chromosome in chrom.sizes

**Problem:** "chr not found in chrom.sizes"
**Solution:** Ensure all chromosomes in bedGraph are in chrom.sizes

## Resources

- [pyBigWig GitHub](https://github.com/deeptools/pyBigWig)
- [UCSC Tools](http://hgdownload.soe.ucsc.edu/admin/exe/)
- [deepTools Documentation](https://deeptools.readthedocs.io/)
- [bigWig Format](https://genome.ucsc.edu/goldenPath/help/bigWig.html)
