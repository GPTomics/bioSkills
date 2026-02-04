# Splicing Quantification - Usage Guide

## Overview
Quantify alternative splicing events from RNA-seq data by calculating PSI (percent spliced in) values. Supports multiple event types including skipped exons, alternative splice sites, mutually exclusive exons, and retained introns.

## Prerequisites
```bash
# SUPPA2 (TPM-based)
pip install suppa

# rMATS-turbo (BAM-based)
pip install rmats-turbo

# Additional dependencies
pip install pandas numpy
```

## Quick Start
Tell your AI agent what you want to do:
- "Quantify exon skipping events from my RNA-seq data"
- "Calculate PSI values for all splicing events in my samples"
- "Generate splicing event quantification from my Salmon TPM files"
- "Run rMATS to quantify splicing from my BAM files"

## Example Prompts

### SUPPA2 Analysis
> "I have transcript TPM values from Salmon. Quantify all alternative splicing events using SUPPA2."

> "Generate PSI values for skipped exon events from my kallisto quantification."

### rMATS Analysis
> "Quantify splicing events from my aligned BAMs using rMATS-turbo."

> "Calculate inclusion levels for mutually exclusive exons across my samples."

### Event-Specific
> "Focus on retained intron events in my RNA-seq data."

> "Quantify alternative 5' and 3' splice site usage."

## What the Agent Will Do
1. Select appropriate tool based on input (TPM vs BAM)
2. Generate splicing event annotations from GTF
3. Calculate PSI values for each event type
4. Filter events by junction read support (minimum 10-20 reads)
5. Output PSI matrices for downstream analysis

## Tips
- SUPPA2 requires transcript-level TPM from Salmon or kallisto
- rMATS-turbo works directly from BAM files
- Minimum 10-20 junction reads recommended for reliable PSI
- PSI values range 0-1 where 1 = fully included, 0 = fully excluded
- Events with PSI consistently near 0 or 1 are constitutive, not alternative

## Related Skills
- differential-splicing - Compare PSI between conditions
- rna-quantification/alignment-free-quant - Generate transcript TPM for SUPPA2
- read-alignment/star-alignment - Align reads with junction detection
