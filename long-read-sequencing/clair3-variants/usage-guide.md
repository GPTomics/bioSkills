# Clair3 Variant Calling Usage Guide

## Overview

Clair3 uses deep learning models trained specifically for long-read sequencing error profiles, achieving state-of-the-art accuracy for germline variant calling from ONT and PacBio data.

## When to Use Clair3

| Scenario | Recommendation |
|----------|---------------|
| ONT germline variants | Clair3 (excellent) |
| PacBio HiFi variants | Clair3 or DeepVariant |
| Somatic variants | ClairS (companion tool) |
| Structural variants | Use Sniffles or cuteSV instead |

## Quick Start Prompts

- "Call variants from my ONT BAM with Clair3"
- "Run Clair3 on PacBio HiFi data"
- "Generate gVCF for joint calling"
- "Filter Clair3 VCF for high-quality SNPs"

## Requirements

```bash
# Install via conda
conda install -c bioconda clair3

# Or via Docker
docker pull hkubal/clair3:latest

# Post-processing
conda install -c bioconda bcftools
```

## Input Requirements

| Parameter | Requirement |
|-----------|-------------|
| BAM | Coordinate-sorted, indexed |
| Reference | Same as used for alignment |
| Coverage | 20-60x recommended |
| Base quality | Higher is better (Q10+) |

## Output Files

| File | Description |
|------|-------------|
| merge_output.vcf.gz | Final variant calls |
| pileup.vcf.gz | Pileup candidates |
| full_alignment_*.vcf.gz | Per-contig calls |
| phased_merge_output.vcf.gz | Phased variants (if enabled) |

## Quality Metrics

| Metric | Good Value |
|--------|------------|
| QUAL | >20 |
| GQ | >30 |
| DP | >10 |

## Related Skills

- **variant-calling/bcftools-basics** - VCF processing
- **long-read-sequencing/alignment** - Input preparation
- **variant-calling/filtering-best-practices** - Quality filtering
