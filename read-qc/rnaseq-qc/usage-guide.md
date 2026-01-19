# RNA-seq QC Usage Guide

This guide covers RNA-seq specific quality control beyond general read quality.

## Key Metrics

### rRNA Contamination
Measure of ribosomal RNA remaining after depletion/selection. Should be <10%.

### Strandedness
Library strand orientation. Must match analysis settings.

### Gene Body Coverage
Even coverage indicates good RNA integrity. 3' bias suggests degradation.

### Transcript Integrity Number (TIN)
Per-transcript measure of degradation. Mean >70 is good.

## Workflow

```
FASTQ
  |
  v
General QC (FastQC)
  |
  v
Alignment
  |
  v
RNA-seq QC <-- This skill
  |
  v
Quantification
```

## Requirements

```bash
# RSeQC
pip install RSeQC

# SortMeRNA
conda install -c bioconda sortmerna

# Picard
conda install -c bioconda picard

# MultiQC
pip install multiqc
```

## Common Issues

### High rRNA
- rRNA depletion failed
- Use SortMeRNA to filter

### Wrong strandedness
- Verify library prep protocol
- Use salmon `-l A` to auto-detect

### 3' bias
- RNA degradation
- Check input RNA quality
- Consider excluding low-TIN samples
