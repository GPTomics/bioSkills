# Tumor Mutational Burden - Usage Guide

## Overview

Calculate tumor mutational burden (TMB) from somatic variant calls with panel-specific normalization to assess immunotherapy eligibility and tumor immunogenicity.

## Prerequisites

```bash
pip install cyvcf2 pandas numpy
```

## Quick Start

Tell your AI agent:
- "Calculate TMB from my somatic VCF"
- "What is the TMB for this tumor using MSK-IMPACT panel size?"
- "Classify TMB as high or low using FDA threshold"
- "Calculate TMB for all samples in this directory"

## Example Prompts

### Basic TMB

> "Calculate TMB from my Mutect2 VCF using FoundationOne CDx panel size"

> "What is the tumor mutational burden for this sample?"

### With Filtering

> "Calculate TMB including only variants with VAF > 5% and depth > 100"

> "Exclude germline variants and calculate TMB"

### Classification

> "Is this sample TMB-high using the FDA 10 mut/Mb threshold?"

> "Compare TMB with MSI status for these samples"

## What the Agent Will Do

1. Load somatic VCF with variant annotations
2. Filter variants (VAF, depth, germline exclusion)
3. Count coding nonsynonymous mutations
4. Divide by panel capture size (megabases)
5. Classify as TMB-High/Low using clinical threshold
6. Optionally compare with MSI status

## Panel Sizes Reference

| Panel | Size (Mb) |
|-------|-----------|
| FoundationOne CDx | 0.8 |
| MSK-IMPACT | 1.14 |
| TruSight Oncology 500 | 1.94 |
| WES (exome) | ~30 |

## Tips

- **Panel size matters**: TMB from different panels is not directly comparable
- **FDA threshold**: 10 mut/Mb for pembrolizumab pan-tumor approval
- **Exclude germline**: Use gnomAD AF > 1% filter
- **VAF filter**: 5% minimum reduces sequencing artifacts
- **MSI-H correlation**: ~80% of MSI-H tumors are also TMB-H
- **Annotation required**: VCF must have consequence annotations (VEP/SnpEff)
