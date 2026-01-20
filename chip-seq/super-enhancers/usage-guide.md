# Super-Enhancer Calling Usage Guide

Identify super-enhancers from H3K27ac ChIP-seq data.

## Prerequisites

```bash
# ROSE
git clone https://github.com/stjude/ROSE.git
# Requires: samtools, R, bedtools

# Alternative: HOMER
conda install -c bioconda homer
```

## Background

Super-enhancers are:
- Large clusters of enhancers (typically > 10kb)
- Marked by H3K27ac, Med1, BRD4
- Control cell identity genes
- Often altered in cancer

## Quick Start

### ROSE (Standard Method)

```bash
python ROSE_main.py \
    -g hg38 \
    -i H3K27ac_peaks.bed \
    -r H3K27ac.bam \
    -o output_directory \
    -s 12500 \  # Stitching distance
    -t 2500     # TSS exclusion zone
```

### Input Requirements

1. **Peak file** (BED or GFF): H3K27ac peaks from MACS3
2. **BAM file**: Aligned H3K27ac ChIP-seq reads
3. **Control BAM** (optional): Input control

### Output Files

- `*_AllEnhancers.table.txt`: All stitched enhancers with signal
- `*_SuperEnhancers.table.txt`: Super-enhancers only
- `*_Enhancers_withSuper.bed`: BED file with SE annotations

## HOMER Alternative

```bash
# Find super-enhancers
findPeaks H3K27ac_tagdir/ -style super \
    -o auto \
    -superSlope 1000 \
    -L 0 \
    -fdr 0.001
```

## Interpretation

### Hockey Stick Plot

ROSE generates a plot showing:
- X-axis: Enhancers ranked by signal
- Y-axis: H3K27ac signal
- Inflection point separates typical enhancers from super-enhancers

### Key Metrics

| Metric | Description |
|--------|-------------|
| Signal | H3K27ac reads in region |
| Rank | Position in ranked list |
| isSuper | 1 if super-enhancer |

## Downstream Analysis

### Assign to Genes

```bash
# Nearest gene assignment
bedtools closest \
    -a SuperEnhancers.bed \
    -b genes.bed \
    > SE_genes.txt
```

### Differential Super-Enhancers

```python
import pandas as pd

# Compare SE signal between conditions
se_a = pd.read_csv('conditionA_SuperEnhancers.txt', sep='\t')
se_b = pd.read_csv('conditionB_SuperEnhancers.txt', sep='\t')

# Find condition-specific SEs
```

## Tips

- Use narrow peaks (not broad) as input
- Increase stitching distance for tissue samples
- Exclude TSS regions to avoid promoter signals
- Validate with independent marks (Med1, BRD4)

## See Also

- [ROSE paper](https://www.cell.com/cell/fulltext/S0092-8674(13)00392-9)
- [Super-enhancer review](https://www.nature.com/articles/nrg.2016.167)
