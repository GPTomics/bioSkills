# Alignment Validation Usage Guide

This guide covers quality control metrics for aligned sequencing data.

## Key Metrics

### Mapping Rate
Percentage of reads successfully aligned. Should be >95% for most experiments.

### Proper Pairing
Percentage of paired reads aligned with correct orientation and insert size. Should be >90%.

### Insert Size
Fragment size distribution. Should match library prep protocol (~300-500bp for WGS).

### Strand Balance
Forward/reverse strand ratio. Should be ~0.5 (balanced).

### GC Bias
Coverage correlation with GC content. Should be flat (no bias).

## Workflow

```
Alignment
    |
    v
Validation <-- This skill
    |
    v
Decision: Pass/Filter/Re-align
    |
    v
Downstream Analysis
```

## Requirements

```bash
# samtools
conda install -c bioconda samtools

# Picard
conda install -c bioconda picard

# pysam
pip install pysam matplotlib
```

## Troubleshooting

### Low mapping rate
- Check reference genome version
- Verify sample species
- Check for contamination

### Poor proper pairing
- Insert size mismatch
- Chimeric reads
- Structural variants

### GC bias
- PCR amplification issues
- Library prep problems
- Consider GC correction
