# ChIP-seq QC Usage Guide

This guide covers quality control metrics for ChIP-seq experiments.

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

## Typical Workflow

1. Align reads and call peaks
2. Calculate FRiP to check enrichment
3. Run phantompeakqualtools for NSC/RSC
4. Check library complexity (NRF, PBC1)
5. Run IDR on replicates
6. Generate fingerprint plots with deepTools

## Tool Requirements

```bash
conda install -c bioconda bedtools samtools phantompeakqualtools idr deeptools
pip install pybedtools pysam
```

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
