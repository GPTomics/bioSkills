# Quality Metrics - Usage Guide

## Overview
Quality metrics assess IMC data quality at acquisition, preprocessing, and segmentation levels to catch issues before downstream analysis.

## Prerequisites
```bash
pip install numpy scipy scikit-image pandas matplotlib
```

## Quick Start
Tell your AI agent what you want to do:
- "Calculate QC metrics for my IMC acquisition"
- "Check signal-to-noise ratio for all channels"
- "Generate a QC report for my experiment"

## Example Prompts

### Signal Quality
> "Calculate signal-to-noise ratio for each channel in my IMC image"

> "Check which channels have low SNR in my data"

### Channel Correlation
> "Calculate pairwise channel correlations to detect spillover"

> "Check for unexpected correlations between my markers"

### Tissue Quality
> "Calculate tissue coverage and fragmentation metrics"

> "Check my image for acquisition artifacts like hot pixels or striping"

### Batch QC
> "Compare QC metrics across all samples in my experiment"

> "Flag outlier samples based on SNR and coverage"

### Comprehensive QC
> "Generate a full QC report for my IMC dataset"

> "Run quality control checks before starting analysis"

## What the Agent Will Do
1. Load preprocessed IMC images
2. Calculate signal-to-noise ratio per channel
3. Compute pairwise channel correlations
4. Assess tissue coverage and fragmentation
5. Detect acquisition artifacts (hot pixels, striping, saturation)
6. Compare metrics across samples for batch effects
7. Generate QC summary report with pass/fail status

## Tips
- Run QC before any downstream analysis
- SNR thresholds: >5 excellent, 3-5 good, 1.5-3 acceptable, <1.5 poor
- Unexpected high correlations (>0.7) may indicate spillover
- Tissue coverage <30% may affect spatial analysis reliability
- Define QC thresholds before looking at data to avoid bias
- Document all QC decisions and excluded samples
- Consider reacquisition for samples failing QC
