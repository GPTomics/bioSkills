# Super-Enhancers - Usage Guide

## Overview

Identify super-enhancers from H3K27ac ChIP-seq data using ROSE or HOMER. Super-enhancers are large clusters of enhancers that control cell identity genes and are often altered in cancer.

## Prerequisites

```bash
# ROSE (standard method)
git clone https://github.com/stjude/ROSE.git
# Requires: samtools, R, bedtools

# HOMER (alternative)
conda install -c bioconda homer
```

## Quick Start

Tell your AI agent what you want to do:
- "Identify super-enhancers from my H3K27ac ChIP-seq data"
- "Run ROSE to find super-enhancers and generate a hockey stick plot"
- "Compare super-enhancers between two conditions"

## Example Prompts

### Basic Super-Enhancer Calling
> "Run ROSE on my H3K27ac peaks to identify super-enhancers"

### Using HOMER
> "Use HOMER findPeaks with -style super to identify super-enhancers"

### Gene Assignment
> "Find which genes are associated with each super-enhancer"

### Differential Analysis
> "Compare super-enhancers between tumor and normal samples"

## What the Agent Will Do

1. Verify H3K27ac peak file and BAM file inputs
2. Run ROSE or HOMER to stitch enhancers and rank by signal
3. Generate hockey stick plot showing inflection point
4. Output super-enhancer BED file and statistics table
5. Optionally assign super-enhancers to nearest genes

## ROSE Workflow

```bash
python ROSE_main.py \
    -g hg38 \
    -i H3K27ac_peaks.bed \
    -r H3K27ac.bam \
    -o output_directory \
    -s 12500 \
    -t 2500
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
findPeaks H3K27ac_tagdir/ -style super \
    -o auto \
    -superSlope 1000 \
    -L 0 \
    -fdr 0.001
```

## Tips

- Use narrow peaks (not broad) as input
- Increase stitching distance for tissue samples
- Exclude TSS regions to avoid promoter signals
- Validate with independent marks (Med1, BRD4)
- Hockey stick plot inflection point separates typical enhancers from super-enhancers

## Downstream Analysis

### Assign to Genes

```bash
bedtools closest \
    -a SuperEnhancers.bed \
    -b genes.bed \
    > SE_genes.txt
```

### Differential Super-Enhancers

```python
import pandas as pd

se_a = pd.read_csv('conditionA_SuperEnhancers.txt', sep='\t')
se_b = pd.read_csv('conditionB_SuperEnhancers.txt', sep='\t')
```

## See Also

- [ROSE paper](https://www.cell.com/cell/fulltext/S0092-8674(13)00392-9)
- [Super-enhancer review](https://www.nature.com/articles/nrg.2016.167)
