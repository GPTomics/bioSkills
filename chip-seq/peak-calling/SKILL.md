---
name: bio-chipseq-peak-calling
description: ChIP-seq peak calling using MACS3 (or MACS2). Call narrow peaks for transcription factors or broad peaks for histone modifications. Supports input control, fragment size modeling, and various output formats including narrowPeak and broadPeak BED files. Use when calling peaks from ChIP-seq alignments.
tool_type: cli
primary_tool: macs3
---

# Peak Calling with MACS3

MACS3 is the actively developed successor to MACS2. Commands are identical except the binary name. MACS2 is in maintenance mode.

## Basic Peak Calling

```bash
# Call peaks with input control (recommended)
macs3 callpeak -t chip.bam -c input.bam -f BAM -g hs -n sample --outdir peaks/

# For MACS2 (legacy), replace 'macs3' with 'macs2' - syntax is identical
```

## Without Input Control

```bash
# Not recommended, but possible
macs3 callpeak -t chip.bam -f BAM -g hs -n sample --outdir peaks/
```

## Narrow Peaks (TF, H3K4me3, H3K27ac)

```bash
macs3 callpeak \
    -t chip.bam \
    -c input.bam \
    -f BAM \
    -g hs \                        # hs=human, mm=mouse, ce=worm, dm=fly
    -n sample_narrow \
    --outdir peaks/ \
    -q 0.05                        # q-value threshold
```

## Broad Peaks (H3K36me3, H3K27me3, H3K9me3)

```bash
macs3 callpeak \
    -t chip.bam \
    -c input.bam \
    -f BAM \
    -g hs \
    -n sample_broad \
    --outdir peaks/ \
    --broad \                      # Broad peak mode
    --broad-cutoff 0.1             # Broad peak q-value
```

## Paired-End Data

```bash
# MACS3 uses BAMPE format for paired-end
macs3 callpeak \
    -t chip.bam \
    -c input.bam \
    -f BAMPE \                     # Paired-end BAM
    -g hs \
    -n sample_pe \
    --outdir peaks/
```

## Multiple Replicates

```bash
# Pool replicates (MACS3 handles internally)
macs3 callpeak \
    -t rep1.bam rep2.bam rep3.bam \
    -c input.bam \
    -f BAM \
    -g hs \
    -n pooled \
    --outdir peaks/
```

## Custom Genome Size

```bash
# For non-model organisms or custom genomes
macs3 callpeak \
    -t chip.bam \
    -c input.bam \
    -f BAM \
    -g 2.7e9 \                     # Effective genome size in bp
    -n sample \
    --outdir peaks/
```

## Common Genome Sizes

| Genome | Flag | Effective Size |
|--------|------|----------------|
| Human | hs | 2.7e9 |
| Mouse | mm | 1.87e9 |
| C. elegans | ce | 9e7 |
| D. melanogaster | dm | 1.2e8 |

## Fixed Fragment Size

```bash
# If modeling fails or for ATAC-seq
macs3 callpeak \
    -t chip.bam \
    -c input.bam \
    -f BAM \
    -g hs \
    --nomodel \                    # Skip model building
    --extsize 200 \                # Fixed extension size
    -n sample \
    --outdir peaks/
```

## Generate Signal Tracks

```bash
# Generate bedGraph and bigWig files
macs3 callpeak \
    -t chip.bam \
    -c input.bam \
    -f BAM \
    -g hs \
    -n sample \
    --outdir peaks/ \
    -B \                           # Generate bedGraph
    --SPMR                         # Signal per million reads

# Convert to bigWig (requires UCSC tools)
sort -k1,1 -k2,2n peaks/sample_treat_pileup.bdg > peaks/sample.sorted.bdg
bedGraphToBigWig peaks/sample.sorted.bdg chrom.sizes peaks/sample.bw
```

## Local Lambda for Broad Marks

```bash
# Recommended for very broad marks
macs3 callpeak \
    -t chip.bam \
    -c input.bam \
    -f BAM \
    -g hs \
    --broad \
    --nolambda \                   # Use local lambda only
    -n sample \
    --outdir peaks/
```

## Cutoff Analysis

```bash
# Test different q-value cutoffs
macs3 callpeak \
    -t chip.bam \
    -c input.bam \
    -f BAM \
    -g hs \
    --cutoff-analysis \            # Generate cutoff analysis file
    -n sample \
    --outdir peaks/
```

## Output Files

| File | Description |
|------|-------------|
| *_peaks.narrowPeak | Peak coordinates (BED6+4) |
| *_peaks.broadPeak | Broad peak coordinates |
| *_summits.bed | Peak summit positions |
| *_model.r | R script for model visualization |
| *_treat_pileup.bdg | Treatment signal (with -B) |
| *_control_lambda.bdg | Control signal (with -B) |

## narrowPeak Format

```
chr1  100  200  peak_1  100  .  5.2  10.5  8.3  50
```
Columns: chr, start, end, name, score, strand, signalValue, pValue, qValue, peak

## Filter Peaks

```bash
# Filter by q-value
awk '$9 > 2' peaks.narrowPeak > peaks.filtered.narrowPeak  # -log10(q) > 2 means q < 0.01

# Sort by signal strength
sort -k7,7nr peaks.narrowPeak > peaks.sorted.narrowPeak
```

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| -t | required | Treatment BAM file(s) |
| -c | none | Control BAM file(s) |
| -f | AUTO | Format (BAM, BAMPE, BED) |
| -g | hs | Genome size |
| -n | NA | Output prefix |
| -q | 0.05 | Q-value cutoff |
| -p | none | P-value cutoff (overrides -q) |
| --broad | false | Broad peak calling |
| --nomodel | false | Skip model building |
| --extsize | 200 | Extension size (with --nomodel) |
| -B | false | Generate bedGraph |
| --SPMR | false | Signal per million reads |

## Related Skills

- peak-annotation - Annotate peaks to genes
- differential-binding - Compare peaks between conditions
- alignment-files - Prepare BAM files
- chipseq-visualization - Visualize peaks
