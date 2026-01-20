# DeepVariant Usage Guide

Deep learning variant calling for high-accuracy germline SNP/indel detection.

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

## Model Types

| Model | Data Type |
|-------|-----------|
| `WGS` | Illumina whole genome |
| `WES` | Illumina exome/targeted |
| `PACBIO` | PacBio HiFi |
| `ONT_R104` | Oxford Nanopore R10.4 |

## Usage Patterns

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
