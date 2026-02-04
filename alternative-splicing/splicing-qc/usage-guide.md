# Splicing QC - Usage Guide

## Overview
Assess RNA-seq data quality specifically for splicing analysis. Evaluates junction saturation, splice site strength, and coverage metrics to determine if data is suitable for reliable splicing quantification.

## Prerequisites
```bash
# RSeQC
pip install rseqc

# MaxEntScan scoring
pip install maxentpy

# Additional dependencies
pip install pysam pandas matplotlib
```

## Quick Start
Tell your AI agent what you want to do:
- "Check if my RNA-seq data has sufficient depth for splicing analysis"
- "Run junction saturation analysis on my BAM files"
- "Evaluate splice site strength for my detected junctions"
- "Assess junction coverage across my samples"

## Example Prompts

### Junction Saturation
> "Run junction saturation analysis to check if I have enough sequencing depth for splicing."

> "Generate saturation curves for junction detection in my samples."

### Junction Quality
> "Classify my detected junctions as known or novel using the annotation."

> "What proportion of my junctions are in the reference annotation?"

### Splice Site Analysis
> "Score the splice site strength for my differential splicing events."

> "Identify weak splice sites that might indicate regulatory splicing."

## What the Agent Will Do
1. Run junction saturation analysis with RSeQC
2. Classify junctions as known/novel using annotation
3. Calculate junction read coverage distribution
4. Score splice site strength with MaxEntScan
5. Generate QC summary with pass/fail recommendations

## Tips
- Junction saturation should reach a plateau for reliable splicing analysis
- If curves are still rising, consider deeper sequencing
- High proportion (>80%) of known junctions indicates good alignment
- 5' splice site scores typically 8-10 bits, 3' sites 8-12 bits
- Weak splice sites (score < 5) may indicate cryptic or regulated splicing
- Note: RSeQC v3.0+ removed the `-s` flag from junction_saturation.py

## Related Skills
- splicing-quantification - Proceed to quantification after QC passes
- read-alignment/star-alignment - Alignment quality affects junction detection
- read-qc/quality-reports - General sequencing QC
