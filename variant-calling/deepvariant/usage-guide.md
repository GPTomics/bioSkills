# DeepVariant - Usage Guide

## Overview

Deep learning variant calling using Google's DeepVariant for high-accuracy germline SNP and indel detection from Illumina, PacBio HiFi, or Oxford Nanopore data.

## Prerequisites

```bash
# Docker (recommended)
docker pull google/deepvariant:1.6.0

# GPU version
docker pull google/deepvariant:1.6.0-gpu

# Singularity
singularity pull docker://google/deepvariant:1.6.0
```

## Quick Start

Tell your AI agent what you want to do:
- "Call variants from my WGS BAM using DeepVariant"
- "Run DeepVariant on exome data with target regions"
- "Set up DeepVariant with GPU acceleration"
- "Generate GVCFs for joint calling with GLnexus"

## Example Prompts

### Basic Variant Calling
> "Call variants from my whole genome BAM using DeepVariant"

> "Run DeepVariant on sample.bam with 16 threads"

### Exome/Targeted
> "Run DeepVariant on my exome data with the capture BED file"

> "Call variants only in my target regions using DeepVariant WES model"

### Long Reads
> "Call variants from my PacBio HiFi reads with DeepVariant"

> "Run DeepVariant with the ONT_R104 model for Nanopore data"

### Multi-Sample
> "Generate GVCFs with DeepVariant for joint calling"

> "Set up DeepVariant + GLnexus for my cohort of 20 samples"

### GPU Acceleration
> "Run DeepVariant with GPU support for faster processing"

## What the Agent Will Do

1. Verify input BAM is aligned, sorted, and indexed
2. Select appropriate model type (WGS, WES, PACBIO, ONT_R104)
3. Set up Docker/Singularity command with volume mounts
4. Run DeepVariant with optimal thread count
5. Generate VCF and optionally GVCF for joint calling
6. Produce variant statistics for quality assessment

## Tips

- Use the correct model type for your data (WGS, WES, PACBIO, ONT_R104)
- Always specify `--regions` for exome/targeted data to save time
- GPU version provides 3-4x speedup for WGS; worthwhile for many samples
- Generate GVCFs with `--output_gvcf` if planning multi-sample joint calling
- Use GLnexus (not GATK) for joint genotyping DeepVariant GVCFs
- Ti/Tv ratio around 2.0-2.1 indicates high-quality WGS SNP calls

## Model Types

| Model | Data Type |
|-------|-----------|
| `WGS` | Illumina whole genome |
| `WES` | Illumina exome/targeted |
| `PACBIO` | PacBio HiFi |
| `ONT_R104` | Oxford Nanopore R10.4 |

## Usage Patterns

### Basic WGS

```bash
docker run -v "${PWD}:/data" google/deepvariant:1.6.0 \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref=/data/reference.fa \
    --reads=/data/sample.bam \
    --output_vcf=/data/output.vcf.gz \
    --output_gvcf=/data/output.g.vcf.gz \
    --num_shards=16
```

### Exome/Targeted Sequencing

```bash
docker run -v "${PWD}:/data" google/deepvariant:1.6.0 \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=WES \
    --ref=/data/reference.fa \
    --reads=/data/exome.bam \
    --regions=/data/targets.bed \
    --output_vcf=/data/output.vcf.gz \
    --num_shards=8
```

### GPU Acceleration

```bash
docker run --gpus all -v "${PWD}:/data" google/deepvariant:1.6.0-gpu \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref=/data/reference.fa \
    --reads=/data/sample.bam \
    --output_vcf=/data/output.vcf.gz
```

### Multi-Sample with GLnexus

1. Generate gVCFs for each sample with `--output_gvcf`
2. Joint genotype with GLnexus:

```bash
docker run -v "${PWD}:/data" quay.io/mlin/glnexus:v1.4.1 \
    /usr/local/bin/glnexus_cli --config DeepVariantWGS \
    /data/*.g.vcf.gz | bcftools view - -Oz -o cohort.vcf.gz
```

## Quality Control

```bash
# Statistics
bcftools stats output.vcf.gz > stats.txt

# Filter low quality
bcftools view -i 'QUAL>20 && FMT/GQ>20' output.vcf.gz -Oz -o filtered.vcf.gz

# Ti/Tv ratio (expect ~2.0-2.1 for WGS)
bcftools stats output.vcf.gz | grep TSTV
```

## Resource Requirements

| Data Type | Memory | Time (30x) |
|-----------|--------|------------|
| WGS | 64 GB | 4-6 hours |
| WES | 32 GB | 30 min |
| WGS + GPU | 32 GB | 1-2 hours |

## See Also

- [DeepVariant GitHub](https://github.com/google/deepvariant)
- [GIAB benchmarking](https://github.com/genome-in-a-bottle)
