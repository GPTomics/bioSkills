# DeepVariant - Usage Guide

## Overview

Deep learning variant calling using Google's DeepVariant for high-accuracy germline SNP and indel detection from Illumina, PacBio HiFi, or Oxford Nanopore data.

## Prerequisites

```bash
# Docker (recommended)
docker pull google/deepvariant:1.6.1

# GPU version
docker pull google/deepvariant:1.6.1-gpu

# DeepTrio (separate image tag)
docker pull google/deepvariant:deeptrio-1.6.1

# Singularity
singularity pull docker://google/deepvariant:1.6.1
```

Input is a sorted, indexed, duplicate-marked BAM/CRAM. Do NOT run BQSR first -- DeepVariant learned its error model from raw base qualities and BQSR slightly lowers its accuracy.

## Quick Start

Tell your AI agent what you want to do:
- "Call variants from my WGS BAM using DeepVariant"
- "Run DeepVariant on exome data with target regions"
- "Set up DeepVariant with GPU acceleration"
- "Generate GVCFs for joint calling with GLnexus"
- "Call de-novo variants in my trio with DeepTrio"
- "Should I use DeepVariant, GATK, or DRAGEN for this cohort?"

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

### Family / Trio
> "Run DeepTrio on my proband-father-mother BAMs and flag de-novo variants"

> "Joint genotype my DeepTrio gVCFs and check Mendelian consistency"

### GPU Acceleration
> "Run DeepVariant with GPU support for faster processing"

## What the Agent Will Do

1. Verify input BAM is aligned, sorted, indexed, and duplicate-marked (and that BQSR was NOT applied)
2. Select the model type matching the sequencing instrument (WGS, WES, PACBIO, ONT_R104, HYBRID_PACBIO_ILLUMINA)
3. Set up the Docker/Singularity command with volume mounts, one-shot or three-stage
4. Run DeepVariant with an appropriate `--num_shards`; use GPU only for call_variants
5. Generate VCF and, when a cohort is planned, gVCF for GLnexus joint calling
6. For families, run DeepTrio and merge the trio gVCFs with GLnexus
7. Produce variant statistics and, when a truth set is available, a stratified hap.py/vcfeval benchmark

## Tips

- Match the model to the instrument, not the analysis goal -- the wrong `--model_type` silently degrades accuracy without erroring
- Do NOT apply GATK hard filters or VQSR to DeepVariant output; the CNN already filtered it (FILTER = PASS / RefCall). Threshold on QUAL/GQ only
- Do NOT run BQSR upstream; it slightly lowers DeepVariant accuracy
- Always specify `--regions` for exome/targeted data to save time
- GPU accelerates only the call_variants stage; make_examples and postprocess_variants are CPU-bound (scale with `--num_shards`)
- Generate gVCFs with `--output_gvcf` if planning multi-sample joint calling
- Use GLnexus (not GATK GenotypeGVCFs) for joint genotyping DeepVariant gVCFs
- Use DeepTrio for de-novo calling; naive trio subtraction produces mostly false de-novos
- Ti/Tv around 2.0-2.1 (WGS) or 3.0-3.3 (WES) indicates high-quality SNP calls
- DeepVariant is trained on GIAB then benchmarked on GIAB; weight held-out/non-GIAB-ancestry performance and validate on characterized material before clinical use

## Model Types

| Model | Data Type |
|-------|-----------|
| `WGS` | Illumina whole genome |
| `WES` | Illumina exome/targeted |
| `PACBIO` | PacBio HiFi |
| `ONT_R104` | Oxford Nanopore R10.4+ |
| `HYBRID_PACBIO_ILLUMINA` | Samples with both HiFi and Illumina |

## Usage Patterns

### Basic WGS

```bash
docker run -v "${PWD}:/data" google/deepvariant:1.6.1 \
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
docker run -v "${PWD}:/data" google/deepvariant:1.6.1 \
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
docker run --gpus all -v "${PWD}:/data" google/deepvariant:1.6.1-gpu \
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

## Related Skills

- variant-calling/variant-calling - Engine-selection table (DeepVariant vs GATK vs DRAGEN) and lightweight bcftools calling
- variant-calling/gatk-variant-calling - GATK HaplotypeCaller alternative with joint calling and VQSR/VETS
- variant-calling/joint-calling - GATK reference-confidence joint genotyping, the alternative to GLnexus
- variant-calling/vcf-basics - View and query VCF output
- variant-calling/filtering-best-practices - Post-call filtering for callers that expose hard-filter annotations (not DeepVariant)
- variant-calling/vcf-statistics - QC metrics for called variants
- long-read-sequencing/clair3-variants - Long-read alternative, especially for ONT R9.4
