# Variant Filtering Best Practices Usage Guide

## Overview

Variant filtering removes false positives while retaining true variants. The approach depends on dataset size, variant type, and analysis goals.

## Quick Start Prompts

- "Apply GATK hard filters to my SNP calls"
- "Filter variants with bcftools for quality and depth"
- "Explain what QD and FS metrics mean"
- "Set up filtering for somatic variant calls"

## Filter Strategy Selection

| Dataset | Recommended Approach |
|---------|---------------------|
| >30 WGS samples | VQSR (machine learning) |
| Small germline | Hard filters |
| Somatic/tumor | Caller-specific + manual |
| Population studies | MAF + HWE + missingness |

## Key Quality Metrics

| Metric | Threshold | Meaning |
|--------|-----------|---------|
| QUAL | >30 | Phred-scaled quality |
| QD | >2 | Quality by depth |
| FS | <60 (SNP) | Strand bias |
| MQ | >40 | Mapping quality |
| DP | 10-500x | Depth (sample-specific) |
| GQ | >20 | Genotype confidence |

## Common Workflow

1. **Separate** SNPs and indels (different thresholds)
2. **Apply** quality filters (QUAL, QD, FS, MQ)
3. **Filter** depth (avoid extremes)
4. **Exclude** problematic regions
5. **Validate** with Ti/Tv ratio

## Requirements

```bash
# bcftools
conda install -c bioconda bcftools

# GATK
conda install -c bioconda gatk4
```

## Related Skills

- **variant-calling/gatk-variant-calling** - VQSR details
- **variant-calling/bcftools-variant-calling** - bcftools usage
