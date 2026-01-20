# ChIP-seq QC - Usage Guide

## Overview

Quality control metrics for ChIP-seq experiments including FRiP, NSC/RSC, IDR, and library complexity measurements to assess enrichment quality and replicate reproducibility.

## Prerequisites

```bash
conda install -c bioconda bedtools samtools phantompeakqualtools idr deeptools
pip install pybedtools pysam
```

## Quick Start

Tell your AI agent what you want to do:
- "Calculate FRiP score for my ChIP-seq sample"
- "Run phantompeakqualtools to get NSC and RSC metrics"
- "Check reproducibility between replicates using IDR"

## Example Prompts

### Enrichment Assessment
> "Calculate the FRiP score for my H3K27ac ChIP-seq experiment"

### Cross-correlation Analysis
> "Run phantompeakqualtools on my BAM file and interpret the NSC/RSC values"

### Replicate Consistency
> "Run IDR analysis on my two ChIP-seq replicates to assess reproducibility"

### Comprehensive QC
> "Generate a full QC report including FRiP, NSC/RSC, and fingerprint plots"

## What the Agent Will Do

1. Calculate FRiP score by counting reads in peaks vs total reads
2. Run phantompeakqualtools to compute NSC and RSC metrics
3. Check library complexity metrics (NRF, PBC1)
4. Run IDR on replicate peak files
5. Generate fingerprint plots with deepTools plotFingerprint
6. Interpret results against ENCODE standards

## Key Metrics

### FRiP (Fraction of Reads in Peaks)
- Proportion of total reads falling within called peaks
- Indicates enrichment strength
- TF ChIP: > 1% minimum, > 5% ideal
- Histone ChIP: > 10% minimum, > 20% ideal

### NSC (Normalized Strand Coefficient)
- Ratio of cross-correlation at fragment length vs background
- Values > 1.1 indicate good enrichment
- Calculated by phantompeakqualtools

### RSC (Relative Strand Coefficient)
- Ratio of fragment-length peak to read-length peak
- Values > 1.0 indicate good enrichment
- More robust than NSC for low-quality samples

### IDR (Irreproducibility Discovery Rate)
- Measures consistency between biological replicates
- Ranks peaks by signal and checks concordance
- > 70% of peaks should be reproducible at IDR < 0.05

## Tips

- Always use input control when calculating enrichment metrics
- Run QC before differential binding analysis to identify problematic samples
- Low FRiP may indicate weak antibody or insufficient sequencing depth
- Poor IDR suggests biological variation or technical issues with one replicate
- Consider re-doing immunoprecipitation if NSC/RSC values are consistently low

## Troubleshooting

### Low FRiP
- Check antibody specificity
- Verify input control quality
- Increase sequencing depth

### Low NSC/RSC
- Fragment size may not match protocol
- Weak enrichment
- Consider re-doing immunoprecipitation

### Poor IDR
- Biological variation between replicates
- Technical issues with one replicate
- Consider calling peaks on pooled data
