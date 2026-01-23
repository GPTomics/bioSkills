# Variant Filtering Best Practices - Usage Guide

## Overview

Variant filtering removes false positives while retaining true variants. The approach depends on dataset size, variant type, and analysis goals.

## Prerequisites

```bash
# bcftools
conda install -c bioconda bcftools

# GATK
conda install -c bioconda gatk4
```

## Quick Start

Tell your AI agent what you want to do:
- "Apply GATK hard filters to my SNP calls"
- "Filter variants with bcftools for quality and depth"
- "Set up VQSR for my whole genome cohort"
- "Filter somatic variants from my tumor-normal calls"

## Example Prompts

### Hard Filtering
> "Apply GATK hard filters to my SNP calls"

> "Filter variants with bcftools for quality and depth"

> "Apply standard hard filter thresholds for indels"

### VQSR
> "Run VQSR on my whole genome cohort with 100 samples"

> "Apply VQSR filtering using HapMap and 1000G resources"

### Understanding Metrics
> "Explain what QD and FS metrics mean in my VCF"

> "Help me interpret the FILTER column annotations"

### Somatic Filtering
> "Set up filtering for somatic variant calls from Mutect2"

> "Filter my tumor-only variant calls to reduce false positives"

### Population Filtering
> "Filter my cohort VCF by MAF and Hardy-Weinberg equilibrium"

> "Remove variants with high missingness across samples"

## What the Agent Will Do

1. Assess dataset characteristics (size, variant type, caller used)
2. Recommend filtering strategy (VQSR vs hard filters)
3. Separate SNPs and indels if using different thresholds
4. Apply appropriate quality filters (QUAL, QD, FS, MQ, DP)
5. Exclude problematic regions if blacklist provided
6. Validate filtering with Ti/Tv ratio and other QC metrics

## Tips

- Always separate SNPs and indels before filtering (different optimal thresholds)
- Ti/Tv ratio of ~2.0-2.1 indicates good quality WGS SNP calls
- Extreme depth values often indicate mapping artifacts; filter both low and high
- VQSR requires many variants to train; use hard filters for exomes or small cohorts
- For somatic variants, use caller-specific filtering recommendations (e.g., Mutect2 FilterMutectCalls)

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

## Related Skills

- **variant-calling/gatk-variant-calling** - VQSR details
- **variant-calling/bcftools-variant-calling** - bcftools usage
