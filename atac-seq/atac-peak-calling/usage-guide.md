# ATAC-seq Peak Calling Usage Guide

## Overview

ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing) uses Tn5 transposase to identify open chromatin regions. Peak calling identifies these accessible regions.

## Key Differences from ChIP-seq

| Aspect | ATAC-seq | ChIP-seq |
|--------|----------|----------|
| Signal source | Tn5 cut sites | Protein binding |
| Control | No input control | Input/IgG required |
| Fragment pattern | Nucleosome periodicity | Smooth enrichment |
| Duplicates | Keep all | Remove PCR duplicates |

## Complete Workflow

```bash
#!/bin/bash
set -euo pipefail

# Variables
BAM=$1
SAMPLE=$(basename $BAM .bam)
GENOME=hs  # hs, mm, dm, ce
OUTDIR=atac_analysis/$SAMPLE

mkdir -p $OUTDIR

# 1. Filter for properly paired, high-quality reads
samtools view -b -f 2 -F 1804 -q 30 $BAM > $OUTDIR/${SAMPLE}.filtered.bam
samtools index $OUTDIR/${SAMPLE}.filtered.bam

# 2. Remove mitochondrial reads
samtools view -h $OUTDIR/${SAMPLE}.filtered.bam | \
    grep -v chrM | \
    samtools view -b > $OUTDIR/${SAMPLE}.noMT.bam
samtools index $OUTDIR/${SAMPLE}.noMT.bam

# 3. Separate NFR and nucleosomal reads
# NFR (<100bp)
samtools view -h $OUTDIR/${SAMPLE}.noMT.bam | \
    awk 'substr($0,1,1)=="@" || ($9>0 && $9<100) || ($9<0 && $9>-100)' | \
    samtools view -b > $OUTDIR/${SAMPLE}.nfr.bam

# Mono-nucleosomal (180-247bp)
samtools view -h $OUTDIR/${SAMPLE}.noMT.bam | \
    awk 'substr($0,1,1)=="@" || ($9>=180 && $9<=247) || ($9<=-180 && $9>=-247)' | \
    samtools view -b > $OUTDIR/${SAMPLE}.mono.bam

# 4. Call peaks on full data
macs3 callpeak \
    -t $OUTDIR/${SAMPLE}.noMT.bam \
    -f BAMPE \
    -g $GENOME \
    -n ${SAMPLE} \
    --outdir $OUTDIR/peaks \
    --nomodel \
    --shift -75 \
    --extsize 150 \
    --keep-dup all \
    -q 0.05 \
    -B \
    --call-summits

# 5. Call peaks on NFR only (for TF binding sites)
macs3 callpeak \
    -t $OUTDIR/${SAMPLE}.nfr.bam \
    -f BAMPE \
    -g $GENOME \
    -n ${SAMPLE}_nfr \
    --outdir $OUTDIR/peaks \
    --nomodel \
    --shift -37 \
    --extsize 75 \
    --keep-dup all \
    -q 0.01 \
    --call-summits

echo "Peaks in $OUTDIR/peaks/"
```

## Understanding Tn5 Offset

Tn5 creates a 9bp duplication at insertion sites:
- Forward read 5' end: +4bp from actual cut
- Reverse read 5' end: -5bp from actual cut

The `--shift` parameter adjusts for this:
```bash
# For paired-end (BAMPE), MACS3 uses fragment centers
# --shift -75 shifts the center to approximate Tn5 cuts

# For single-end, shift read 5' ends
--shift -4  # Forward reads
--shift +5  # Reverse reads (handled by extension)
```

## Fragment Size Selection

ATAC-seq produces characteristic fragment sizes:

| Fragment Size | Origin | Analysis Use |
|---------------|--------|--------------|
| <100 bp | Nucleosome-free | TF binding, footprinting |
| 180-247 bp | Mono-nucleosome | Nucleosome positioning |
| 315-473 bp | Di-nucleosome | Chromatin structure |
| 558-615 bp | Tri-nucleosome | Chromatin structure |

```bash
# Extract specific sizes with samtools
# NFR
samtools view -h in.bam | awk 'BEGIN{OFS="\t"} /^@/ || ($9>0 && $9<100) || ($9<0 && $9>-100)' | samtools view -b > nfr.bam

# Mono-nucleosomal
samtools view -h in.bam | awk 'BEGIN{OFS="\t"} /^@/ || ($9>=180 && $9<=247) || ($9<=-180 && $9>=-247)' | samtools view -b > mono.bam
```

## Quality-Based Peak Filtering

```bash
# Filter by q-value
awk '$9 >= 2' peaks.narrowPeak > peaks.q01.bed  # -log10(0.01) = 2

# Filter by signal strength
awk '$7 >= 10' peaks.narrowPeak > peaks.high_signal.bed

# Size filter
awk '$3-$2 >= 100 && $3-$2 <= 10000' peaks.narrowPeak > peaks.size_filtered.bed
```

## Blacklist Filtering

Remove problematic regions:

```bash
# Download ENCODE blacklist
wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
gunzip hg38-blacklist.v2.bed.gz

# Filter peaks
bedtools intersect -v -a peaks.narrowPeak -b hg38-blacklist.v2.bed > peaks.filtered.bed
```

## Consensus Peaks Across Samples

```bash
# Merge overlapping peaks
cat sample1_peaks.narrowPeak sample2_peaks.narrowPeak sample3_peaks.narrowPeak | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i - > consensus_peaks.bed

# Or use bedtools multiinter for more control
bedtools multiinter -i sample1.narrowPeak sample2.narrowPeak sample3.narrowPeak | \
    awk '$4 >= 2' | \  # Present in at least 2 samples
    bedtools merge -i - > consensus_peaks.bed
```

## Genrich Alternative

Genrich is designed specifically for ATAC-seq:

```bash
# Genrich ATAC-seq mode
Genrich -t sample.bam \
    -o peaks.narrowPeak \
    -j \                    # ATAC-seq mode
    -r \                    # Remove PCR duplicates
    -e chrM \               # Exclude mitochondria
    -q 0.05
```

## Peak Statistics

```python
import pandas as pd

peaks = pd.read_csv('peaks.narrowPeak', sep='\t', header=None,
    names=['chr', 'start', 'end', 'name', 'score', 'strand', 'signal', 'pval', 'qval', 'summit'])

print(f"Total peaks: {len(peaks)}")
print(f"Median peak width: {(peaks['end'] - peaks['start']).median():.0f}")
print(f"Total accessible bp: {(peaks['end'] - peaks['start']).sum():,}")
print(f"Mean signal: {peaks['signal'].mean():.1f}")
```

## Troubleshooting

### Too Few Peaks
- Check alignment quality and mapping rate
- Verify fragment size distribution shows nucleosome pattern
- Try lower q-value threshold

### Too Many Peaks
- Increase q-value threshold
- Check for high mitochondrial content
- Verify blacklist filtering

### Noisy Signal
- Filter for higher MAPQ
- Remove duplicates more stringently
- Check for contamination
