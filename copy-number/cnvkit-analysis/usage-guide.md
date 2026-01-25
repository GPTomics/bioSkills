# CNVkit Analysis Usage Guide

## Overview

CNVkit is the standard tool for detecting copy number variants from targeted sequencing (exome, gene panels). It uses both on-target and off-target reads to infer copy number across the genome.

## Prerequisites

```bash
conda install -c bioconda cnvkit
# or
pip install cnvkit
```

Dependencies: R with DNAcopy package (for CBS segmentation).

## Quick Start

Tell your AI agent what you want to do:
- "Run CNVkit on my tumor-normal exome pair to call copy number variants"
- "Build a panel of normals from my control samples for CNV calling"
- "Export my CNVkit segments to VCF format for downstream analysis"
- "Check quality metrics for my CNVkit run and identify noisy samples"

## Input Requirements

1. **BAM files**: Aligned, sorted, indexed
2. **Target BED**: Capture regions (exome targets)
3. **Reference FASTA**: Same as used for alignment
4. **Gene annotations** (optional): refFlat.txt for gene names

## Workflow Selection

| Scenario | Command |
|----------|---------|
| Tumor-normal pair | `batch tumor.bam --normal normal.bam ...` |
| Multiple normals → reference | `batch --normal *.bam ...` then `batch tumor.bam --reference ref.cnn` |
| Tumor-only (flat reference) | `batch tumor.bam --targets ... --fasta ...` |
| WGS | Add `--method wgs` |
| Germline | Use `--normal` samples only |

## Complete Tumor-Normal Pipeline

```bash
# Variables
TUMOR=tumor.bam
NORMAL=normal.bam
TARGETS=exome_targets.bed
FASTA=reference.fa
REFFLAT=refFlat.txt
OUTDIR=cnvkit_results

mkdir -p $OUTDIR

# Run batch command
cnvkit.py batch $TUMOR \
    --normal $NORMAL \
    --targets $TARGETS \
    --annotate $REFFLAT \
    --fasta $FASTA \
    --output-reference $OUTDIR/reference.cnn \
    --output-dir $OUTDIR

# Call with purity adjustment (if known)
cnvkit.py call $OUTDIR/tumor.cns \
    --purity 0.7 \
    -o $OUTDIR/tumor.call.cns

# Generate plots
cnvkit.py scatter $OUTDIR/tumor.cnr -s $OUTDIR/tumor.cns -o $OUTDIR/tumor_scatter.png
cnvkit.py diagram $OUTDIR/tumor.cnr -s $OUTDIR/tumor.cns -o $OUTDIR/tumor_diagram.pdf
```

## Building a Panel of Normals

For best results, use 5-10+ matched normal samples:

```bash
# Step 1: Generate reference from normals
cnvkit.py batch \
    --normal normal1.bam normal2.bam normal3.bam normal4.bam normal5.bam \
    --targets targets.bed \
    --annotate refFlat.txt \
    --fasta reference.fa \
    --output-reference panel_of_normals.cnn

# Step 2: Process tumors with PON
for tumor in tumor*.bam; do
    sample=$(basename $tumor .bam)
    cnvkit.py batch $tumor \
        --reference panel_of_normals.cnn \
        --output-dir results/$sample
done
```

## Interpreting Results

### CNR file (copy ratios)
- Each row is a bin (target or off-target region)
- `log2` column: log2 ratio vs reference (0 = diploid)
- `weight`: confidence weight

### CNS file (segments)
- Merged adjacent bins with similar log2
- `cn`: absolute copy number (after calling)
- `probes`: number of bins in segment

### Thresholds for calling:
| log2 | CN State | Interpretation |
|------|----------|----------------|
| < -1.1 | 0 | Homozygous deletion |
| -1.1 to -0.25 | 1 | Heterozygous deletion |
| -0.25 to 0.2 | 2 | Diploid |
| 0.2 to 0.7 | 3 | Single copy gain |
| > 0.7 | 4+ | Amplification |

## Quality Control

```bash
# Check metrics for all samples
cnvkit.py metrics *.cnr -s *.cns

# Good quality indicators:
# - Median absolute deviation (MAD) < 0.5
# - Biweight midvariance < 0.5
# - Few bins with extreme values

# Check inferred sex
cnvkit.py sex *.cnr *.cnn
```

## Common Issues

### Noisy data
```bash
# Increase smoothing during segmentation
cnvkit.py segment sample.cnr --smooth-cbs -o sample.cns

# Or use HMM which is more robust
cnvkit.py segment sample.cnr --method hmm -o sample.cns
```

### Low coverage samples
```bash
# Lower minimum coverage thresholds
cnvkit.py batch tumor.bam --normal normal.bam \
    --targets targets.bed --fasta reference.fa \
    --target-avg-size 200
```

### Wrong ploidy/purity
```bash
# Estimate from data
cnvkit.py call sample.cns --center median -o sample.call.cns

# Or specify if known
cnvkit.py call sample.cns --purity 0.7 --ploidy 2 -o sample.call.cns
```

## Integration with Other Tools

### Export for downstream analysis
```bash
# For GISTIC2 (identifying recurrent CNVs)
cnvkit.py export seg *.cns -o cohort.seg

# For integration with SVs
cnvkit.py export vcf sample.cns -o sample.cnv.vcf
```

### Python analysis
```python
import cnvlib
import pandas as pd

# Load and filter
cns = cnvlib.read('sample.cns')
amplified = cns[(cns['log2'] > 0.5) & (cns['cn'] >= 4)]
deleted = cns[(cns['log2'] < -0.5) & (cns['cn'] <= 1)]

# Genes in amplified regions
print(amplified[['chromosome', 'start', 'end', 'gene', 'log2', 'cn']])
```

## Example Prompts

> "Run CNVkit on my tumor-normal exome pair and call copy number variants"

> "Build a panel of normals from my 10 control samples for CNV calling"

> "Export my CNVkit segments to VCF format for integration with SNV calls"

> "Check the quality metrics for my CNVkit run and identify noisy samples"
