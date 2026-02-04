# ctDNA Mutation Detection - Usage Guide

## Overview
Detect somatic mutations in circulating tumor DNA at low variant allele fractions. Uses specialized callers optimized for cfDNA with UMI-based error suppression.

## Prerequisites
```bash
# VarDict
conda install -c bioconda vardict-java

# UMI-VarCal (optional, best specificity)
pip install umi-varcal

# Dependencies
pip install pysam pandas vcfpy
```

## Quick Start
Tell your AI agent what you want to do:
- "Detect mutations at 0.5% VAF in my cfDNA panel data"
- "Call variants from my UMI-consensus BAM"
- "Filter out CHIP mutations from my ctDNA calls"
- "Track specific mutations across serial samples"

## Example Prompts

### Variant Calling
> "Run VarDict to detect mutations at 0.5% VAF from my targeted panel."

> "Call variants from my UMI-consensus BAM with high sensitivity."

### CHIP Filtering
> "Filter out potential CHIP variants (DNMT3A, TET2, ASXL1, etc.)."

> "Separate likely somatic mutations from age-related clonal hematopoiesis."

### Mutation Tracking
> "Track these known mutations across my serial plasma samples."

## What the Agent Will Do
1. Call variants using cfDNA-optimized caller
2. Filter by VAF and read support thresholds
3. Annotate variants with clinical information
4. Flag potential CHIP variants
5. Track specific mutations if requested

## Tips
- Use targeted panels or WES for mutation detection (NOT sWGS)
- Reliable detection above 0.5% VAF with UMIs
- VarDict is highly sensitive but may have more false positives
- UMI-VarCal offers best specificity with UMI data
- Always filter for CHIP (DNMT3A, TET2, ASXL1 common in elderly)
- smCounter2 has 0.5% VAF detection limit

## Related Skills
- cfdna-preprocessing - Preprocess with UMI consensus
- tumor-fraction-estimation - Estimate overall tumor burden
- longitudinal-monitoring - Track mutations over time
- variant-calling/variant-calling - General variant calling concepts
