# Clair3 Variant Calling - Usage Guide

## Overview
Clair3 uses deep learning models trained specifically for long-read sequencing error profiles, achieving state-of-the-art accuracy for germline variant calling from ONT and PacBio data.

## Prerequisites
```bash
# Install via conda
conda install -c bioconda clair3 bcftools

# Or via Docker
docker pull hkubal/clair3:latest
```

## Quick Start
Tell your AI agent what you want to do:
- "Call germline variants from my ONT BAM using Clair3"
- "Run Clair3 on PacBio HiFi data"

## Example Prompts

### Basic Variant Calling
> "Call variants from aligned.bam using Clair3 with the ONT R10.4.1 model"

> "Run Clair3 on my PacBio HiFi alignment and output to variants.vcf"

### With gVCF for Joint Calling
> "Generate gVCF output from Clair3 for later joint calling with other samples"

### Filtering
> "Call variants with Clair3 and filter for high-quality SNPs with QUAL >= 20"

### Specific Regions
> "Run Clair3 on chromosome 21 only for a quick test"

### Phasing
> "Call variants with Clair3 and enable phasing"

## What the Agent Will Do
1. Verify BAM is sorted and indexed
2. Select appropriate Clair3 model for your platform (ONT/HiFi)
3. Run variant calling with suitable parameters
4. Apply quality filtering as requested
5. Index output VCF

## When to Use Clair3

| Scenario | Recommendation |
|----------|---------------|
| ONT germline variants | Clair3 (excellent) |
| PacBio HiFi variants | Clair3 or DeepVariant |
| Somatic variants | ClairS (companion tool) |
| Structural variants | Use Sniffles or cuteSV instead |

## Input Requirements

| Parameter | Requirement |
|-----------|-------------|
| BAM | Coordinate-sorted, indexed |
| Reference | Same as used for alignment |
| Coverage | 20-60x recommended |
| Base quality | Higher is better (Q10+) |

## Quality Metrics

| Metric | Good Value |
|--------|------------|
| QUAL | >20 |
| GQ | >30 |
| DP | >10 |

## Tips
- Use matching platform model (ONT vs HiFi) - wrong model significantly impacts accuracy
- 20-60x coverage gives best results; below 20x reduces sensitivity
- Enable `--gvcf` for joint calling across multiple samples
- For somatic variants, use ClairS instead of Clair3
- Docker is often easier than conda installation

## Resources
- [Clair3 GitHub](https://github.com/HKU-BAL/Clair3)
- [Clair3 Models](https://github.com/nanoporetech/rerio)
