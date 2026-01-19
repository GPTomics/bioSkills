# TF Footprinting Usage Guide

## Overview

Transcription factor footprinting leverages the protection from Tn5 cutting provided by DNA-bound proteins. When a TF binds DNA, it creates a characteristic "footprint" - a dip in ATAC-seq signal at the binding site flanked by high signal where Tn5 can cut.

## Biological Principle

```
High signal    Low signal    High signal
(accessible)   (protected)   (accessible)
   ____          ____          ____
  |    |        |    |        |    |
  |    |________|    |________|    |
                 TF
              binding
```

## Requirements

Footprinting requires:
- **High depth**: >50M uniquely mapped reads
- **NFR reads**: Filter for fragments <100bp
- **Peak regions**: Accessible chromatin regions
- **TF motifs**: JASPAR, HOCOMOCO, or custom

## Complete TOBIAS Pipeline

```bash
#!/bin/bash
set -euo pipefail

BAM=$1
PEAKS=$2
GENOME=$3
MOTIFS=$4
SAMPLE=$(basename $BAM .bam)
OUTDIR=footprinting/$SAMPLE

mkdir -p $OUTDIR

# 1. Filter BAM for NFR
echo "Filtering for NFR reads..."
samtools view -h $BAM | \
    awk 'substr($0,1,1)=="@" || ($9>0 && $9<100) || ($9<0 && $9>-100)' | \
    samtools view -b > $OUTDIR/${SAMPLE}_nfr.bam
samtools index $OUTDIR/${SAMPLE}_nfr.bam

# 2. Correct Tn5 bias
echo "Correcting Tn5 bias..."
tobias ATACorrect \
    --bam $OUTDIR/${SAMPLE}_nfr.bam \
    --genome $GENOME \
    --peaks $PEAKS \
    --outdir $OUTDIR/corrected \
    --cores 8

# 3. Calculate footprint scores
echo "Calculating footprint scores..."
tobias FootprintScores \
    --signal $OUTDIR/corrected/${SAMPLE}_nfr_corrected.bw \
    --regions $PEAKS \
    --output $OUTDIR/${SAMPLE}_footprints.bw \
    --cores 8

# 4. Detect TF binding
echo "Detecting TF binding..."
tobias BINDetect \
    --motifs $MOTIFS \
    --signals $OUTDIR/${SAMPLE}_footprints.bw \
    --genome $GENOME \
    --peaks $PEAKS \
    --outdir $OUTDIR/bindetect \
    --cores 8

echo "Done! Results in $OUTDIR/bindetect/"
```

## Differential Footprinting

Compare TF activity between conditions:

```bash
# Process both conditions first
for cond in condition1 condition2; do
    tobias ATACorrect --bam ${cond}.bam --genome genome.fa --peaks peaks.bed \
        --outdir ${cond}_corrected --cores 8
    tobias FootprintScores --signal ${cond}_corrected/*_corrected.bw \
        --regions peaks.bed --output ${cond}_footprints.bw --cores 8
done

# Differential analysis
tobias BINDetect \
    --motifs JASPAR_motifs.pfm \
    --signals condition1_footprints.bw condition2_footprints.bw \
    --genome genome.fa \
    --peaks peaks.bed \
    --outdir differential \
    --cond_names condition1 condition2 \
    --cores 8
```

## Interpreting TOBIAS Output

### bindetect_results.txt
Main output with per-TF statistics:

| Column | Description |
|--------|-------------|
| name | TF name |
| motif_id | JASPAR ID |
| n_detected | Number of binding sites |
| mean_score | Average footprint score |
| differential_score | Difference between conditions |
| pvalue | Statistical significance |

### Per-TF BED files
- `*_bound.bed`: Predicted bound sites
- `*_unbound.bed`: Predicted unbound sites

## Visualizing Footprints

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pyBigWig

def plot_tf_footprint(footprint_bw, bound_bed, tf_name, output):
    '''Plot aggregate footprint for a TF.'''
    bw = pyBigWig.open(footprint_bw)

    signals = []
    for line in open(bound_bed):
        fields = line.strip().split('\t')
        chrom, start, end = fields[0], int(fields[1]), int(fields[2])
        center = (start + end) // 2

        try:
            vals = bw.values(chrom, center - 100, center + 100)
            if vals:
                signals.append(vals)
        except:
            continue

    avg_signal = np.nanmean(signals, axis=0)
    x = np.arange(-100, 100)

    plt.figure(figsize=(8, 4))
    plt.fill_between(x, avg_signal, alpha=0.5)
    plt.plot(x, avg_signal, 'b-', linewidth=2)
    plt.axvline(0, color='red', linestyle='--', alpha=0.5)
    plt.xlabel('Distance from motif center (bp)')
    plt.ylabel('Footprint score')
    plt.title(f'{tf_name} Aggregate Footprint (n={len(signals)})')
    plt.tight_layout()
    plt.savefig(output, dpi=150)
    plt.close()

# Usage
plot_tf_footprint('footprints.bw', 'CTCF_bound.bed', 'CTCF', 'ctcf_footprint.png')
```

## Quality Assessment

### Good Footprint Indicators
- Clear V-shaped dip at motif center
- Symmetric shoulders
- Many bound sites
- High footprint scores

### Poor Footprint Indicators
- Flat or noisy signal
- Asymmetric pattern
- Few detected sites

```python
def assess_footprint_quality(signal):
    '''Score footprint quality.'''
    center = len(signal) // 2
    shoulder_height = np.mean([signal[:20].mean(), signal[-20:].mean()])
    center_depth = signal[center-5:center+5].mean()

    depth = shoulder_height - center_depth
    symmetry = 1 - abs(signal[:center].mean() - signal[center:].mean()) / shoulder_height

    return {
        'depth': depth,
        'symmetry': symmetry,
        'quality': 'good' if depth > 0.5 and symmetry > 0.8 else 'poor'
    }
```

## Alternative Tools

### HINT-ATAC (RGT)
```bash
rgt-hint footprinting --atac-seq --organism hg38 \
    --output-prefix sample sample.bam peaks.bed

rgt-hint differential --organism hg38 --output-location diff/ \
    --mpbs-file motifs.bed cond1.bam cond2.bam
```

### Wellington (pyDNase)
```python
from pyDNase import GenomicIntervalSet
# Designed for DNase-seq but works for ATAC
```

## Troubleshooting

### No Footprints Detected
- Check read depth (need >50M)
- Verify NFR filtering
- Check peak quality

### All TFs Show Same Pattern
- May indicate poor Tn5 bias correction
- Check input motif quality
- Verify correct genome

### Weak Footprints
- Consider merging replicates
- Use only strongest peaks
- Check for high background
