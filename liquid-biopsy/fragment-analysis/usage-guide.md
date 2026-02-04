# Fragment Analysis - Usage Guide

## Overview
Analyze cfDNA fragment size distributions and fragmentomics patterns for cancer detection. Extract nucleosome positioning signals and DELFI-style fragmentation profiles.

## Prerequisites
```bash
# FinaleToolkit
pip install finaletoolkit

# Griffin (optional, for nucleosome profiling)
pip install griffin

# Dependencies
pip install pysam numpy pandas matplotlib
```

## Quick Start
Tell your AI agent what you want to do:
- "Analyze fragment size distribution for tumor signal"
- "Calculate short-to-long fragment ratios across the genome"
- "Run DELFI-style fragmentomics analysis"
- "Profile nucleosome positioning from my cfDNA"

## Example Prompts

### Fragment Size Analysis
> "Calculate the short (100-150bp) to long (151-220bp) fragment ratio for my sample."

> "Generate genome-wide fragmentation profiles in 5Mb bins."

### FinaleToolkit
> "Run FinaleToolkit to replicate DELFI-style analysis."

> "Calculate GC-corrected fragmentation ratios."

### Nucleosome Profiling
> "Analyze nucleosome accessibility around transcription start sites."

## What the Agent Will Do
1. Extract fragment sizes from BAM files
2. Calculate short/long fragment ratios
3. Generate genome-wide fragmentation profiles
4. Apply GC correction if requested
5. Compare patterns to healthy reference

## Tips
- DELFI is a commercial company, NOT software - use FinaleToolkit (MIT license)
- FinaleToolkit 0.7.1+ is 50x faster than original DELFI approach
- Short fragments (100-150bp) are enriched in ctDNA
- Normal cfDNA peaks at ~167bp (mononucleosome)
- Griffin 0.2.0 is useful for tissue deconvolution

## Related Skills
- cfdna-preprocessing - Preprocess before fragment analysis
- tumor-fraction-estimation - Complement with CNV-based estimation
- methylation-based-detection - Alternative detection approach
