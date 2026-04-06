---
name: bio-variant-calling-deepvariant
description: Deep learning-based variant calling with Google DeepVariant. Provides high accuracy for germline SNPs and indels from Illumina, PacBio, and ONT data. Use when calling variants with DeepVariant deep learning caller or when highest germline calling accuracy is required.
tool_type: cli
primary_tool: DeepVariant
---

## Version Compatibility

Reference examples tested with: DeepVariant 1.6+, GLnexus 1.4+, bcftools 1.19+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `docker run google/deepvariant:<tag> --version` to confirm build
- `bcftools --version` and `bcftools --help` to confirm flags

If code throws errors, introspect the installed container and adapt the example
to match the actual API rather than retrying.

# DeepVariant Variant Calling

## How DeepVariant Works

DeepVariant reframes variant calling as an image classification problem rather than a
statistical genotyping problem. For each candidate variant site, the make_examples step
encodes the local read pileup as a 100x221x6-channel image tensor. The six channels
encode: read base identity, base quality, mapping quality, strand orientation, read
support for the variant allele, and reference base mismatch. A convolutional neural
network (CNN) then classifies each pileup image into one of three genotype classes:
homozygous reference, heterozygous, or homozygous alternate.

This image-based approach is why DeepVariant outperforms statistical callers in
difficult genomic contexts such as homopolymer runs, tandem repeats, and low-complexity
regions -- the CNN learns visual patterns in pileup geometry that heuristic filters
miss. Separate models are trained for each sequencing platform because read error
profiles differ substantially (e.g., Illumina substitution errors vs. ONT indel errors
in homopolymers).

DeepVariant calls germline variants only. For somatic variant calling, use DeepSomatic,
a separate tool from the same team.

## Installation

```bash
docker pull google/deepvariant:1.6.1

# GPU support (NVIDIA GPU + nvidia-container-toolkit required)
docker pull google/deepvariant:1.6.1-gpu
```

Singularity alternative:

```bash
singularity pull docker://google/deepvariant:1.6.1
```

## Model Selection Guide

| Model | Data Type | Notes |
|-------|-----------|-------|
| `WGS` | Illumina short-read WGS | Default model; trained on 30-50x PCR-free data |
| `WES` | Illumina exome/targeted | Must supply `--regions` BED for efficiency; without it, wastes time scanning non-target regions |
| `PACBIO` | PacBio HiFi only | Not trained on CLR reads; CLR error profile is fundamentally different |
| `ONT_R104` | ONT R10.4+ chemistry | Accuracy lower than HiFi model; R9.4 reads perform poorly with this model |
| `HYBRID_PACBIO_ILLUMINA` | Mixed platforms | Emerging use case for samples with both HiFi and Illumina data |

Model selection has a large effect on accuracy. Using the wrong model (e.g., WGS model
on HiFi data) silently degrades results because the CNN expects platform-specific error
patterns in the pileup images.

## When to Use DeepVariant

**Best choice when:**
- Highest germline accuracy is the primary goal
- GPU resources are available (or CPU time is acceptable)
- Single-sample calling or small cohort with GLnexus for joint genotyping
- Data is from a supported platform (Illumina, HiFi, ONT R10.4+)

**Consider GATK HaplotypeCaller instead when:**
- Joint calling across large cohorts (GATK GenomicsDB + GenotypeGVCFs scales better)
- VQSR-based filtering is needed (DeepVariant QUAL scores are CNN confidence, not amenable to VQSR)
- Clinical pipeline requires established GATK validation and regulatory precedent

**Consider Clair3 instead when:**
- Long-read data where speed matters more than marginal accuracy gains
- ONT data specifically (Clair3 has strong ONT-specific models)
- Resource-constrained environments without GPU access

## Basic Usage

```bash
docker run -v "${PWD}:/input" -v "${PWD}/output:/output" \
    google/deepvariant:1.6.1 \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref=/input/reference.fa \
    --reads=/input/sample.bam \
    --output_vcf=/output/sample.vcf.gz \
    --output_gvcf=/output/sample.g.vcf.gz \
    --num_shards=16
```

Always generate gVCFs (`--output_gvcf`) even for single samples -- they enable
downstream joint calling without re-running DeepVariant.

## Step-by-Step Workflow

For more control over intermediate outputs, run each stage separately:

### Step 1: Make Examples

```bash
docker run -v "${PWD}:/data" google/deepvariant:1.6.1 \
    /opt/deepvariant/bin/make_examples \
    --mode calling \
    --ref /data/reference.fa \
    --reads /data/sample.bam \
    --examples /data/examples.tfrecord.gz \
    --gvcf /data/gvcf.tfrecord.gz
```

### Step 2: Call Variants

```bash
docker run -v "${PWD}:/data" google/deepvariant:1.6.1 \
    /opt/deepvariant/bin/call_variants \
    --outfile /data/call_variants.tfrecord.gz \
    --examples /data/examples.tfrecord.gz \
    --checkpoint /opt/models/wgs/model.ckpt
```

### Step 3: Postprocess Variants

```bash
docker run -v "${PWD}:/data" google/deepvariant:1.6.1 \
    /opt/deepvariant/bin/postprocess_variants \
    --ref /data/reference.fa \
    --infile /data/call_variants.tfrecord.gz \
    --outfile /data/output.vcf.gz \
    --gvcf_outfile /data/output.g.vcf.gz \
    --nonvariant_site_tfrecord_path /data/gvcf.tfrecord.gz
```

## GPU Acceleration

GPU acceleration primarily benefits the call_variants step (CNN inference). The
make_examples and postprocess_variants steps are CPU-bound and benefit more from
`--num_shards` parallelism.

```bash
docker run --gpus all -v "${PWD}:/data" \
    google/deepvariant:1.6.1-gpu \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref=/data/reference.fa \
    --reads=/data/sample.bam \
    --output_vcf=/data/output.vcf.gz \
    --num_shards=16
```

## PacBio HiFi Calling

```bash
docker run -v "${PWD}:/data" google/deepvariant:1.6.1 \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=PACBIO \
    --ref=/data/reference.fa \
    --reads=/data/hifi_aligned.bam \
    --output_vcf=/data/hifi_variants.vcf.gz \
    --num_shards=16
```

HiFi reads achieve Q30+ per-read accuracy, giving DeepVariant cleaner pileup images.
The PACBIO model is not suitable for CLR reads (Q10-Q15 accuracy).

## ONT Calling

```bash
docker run -v "${PWD}:/data" google/deepvariant:1.6.1 \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=ONT_R104 \
    --ref=/data/reference.fa \
    --reads=/data/ont_aligned.bam \
    --output_vcf=/data/ont_variants.vcf.gz \
    --num_shards=16
```

R10.4+ chemistry substantially reduces systematic indel errors in homopolymers compared
to R9.4. For R9.4 data, consider Clair3 which has dedicated R9.4 models.

## Exome/Targeted Sequencing

```bash
docker run -v "${PWD}:/data" google/deepvariant:1.6.1 \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=WES \
    --ref=/data/reference.fa \
    --reads=/data/exome.bam \
    --regions=/data/targets.bed \
    --output_vcf=/data/exome_variants.vcf.gz \
    --num_shards=8
```

The `--regions` flag is not strictly required but omitting it causes DeepVariant to
scan the entire genome, wasting hours on off-target reads with minimal coverage.

## Joint Calling with GLnexus

For multi-sample cohorts, use per-sample gVCFs merged with GLnexus:

```bash
for bam in *.bam; do
    sample=$(basename $bam .bam)
    docker run -v "${PWD}:/data" google/deepvariant:1.6.1 \
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=WGS \
        --ref=/data/reference.fa \
        --reads=/data/$bam \
        --output_vcf=/data/${sample}.vcf.gz \
        --output_gvcf=/data/${sample}.g.vcf.gz \
        --num_shards=16
done

docker run -v "${PWD}:/data" quay.io/mlin/glnexus:v1.4.1 \
    /usr/local/bin/glnexus_cli \
    --config DeepVariantWGS \
    /data/*.g.vcf.gz \
    | bcftools view - -Oz -o cohort.vcf.gz
```

### GLnexus Configuration

| Config | Use Case | Notes |
|--------|----------|-------|
| `DeepVariantWGS` | Illumina WGS gVCFs | Default for most WGS cohorts |
| `DeepVariantWES` | Illumina exome gVCFs | Tuned for higher-depth, narrower-region calling |
| `DeepVariant_unfiltered` | Keep all variant sites | Useful for research exploration; produces more false positives |

GLnexus performance is driven by the call confidence distribution across sites, not
cohort size per se. A cohort of 100 samples with clean 30x WGS merges faster than 20
samples with noisy 10x data. GLnexus scales well to thousands of samples.

## Output Quality Metrics

```bash
bcftools stats output.vcf.gz > stats.txt

# Ti/Tv ratio: expect ~2.0-2.1 for WGS, ~2.8-3.3 for WES
bcftools stats output.vcf.gz | grep TSTV

# Filter by quality (QUAL is CNN confidence, GQ is genotype quality)
bcftools view -i 'QUAL>20 && FMT/GQ>20' output.vcf.gz -Oz -o filtered.vcf.gz
```

## Benchmarking Against Truth Set

```bash
docker run -v "${PWD}:/data" jmcdani20/hap.py:latest \
    /opt/hap.py/bin/hap.py \
    /data/HG002_GRCh38_truth.vcf.gz \
    /data/deepvariant_output.vcf.gz \
    -r /data/reference.fa \
    -o /data/benchmark \
    --threads 16
```

## Comparison with Other Callers

Benchmark numbers below are approximate, derived from GIAB HG002/HG003/HG004 truth
sets on GRCh38. Exact values vary by sample, coverage, and version.

| Caller | SNP F1 | Indel F1 | Speed (30x WGS) | Notes |
|--------|--------|----------|------------------|-------|
| DeepVariant | ~0.999 | ~0.993 | ~4-6 hrs CPU, ~1-2 hrs GPU | Highest accuracy; slow without GPU |
| GATK HC | ~0.999 | ~0.989 | ~4-8 hrs CPU | Strong ecosystem; joint calling pipeline |
| Strelka2 | ~0.998 | ~0.960 | ~1-2 hrs CPU | Fastest; no longer actively maintained |
| Clair3 | ~0.998 | ~0.980 | ~8 hrs (50x ONT) | Strong for long reads; active development |

Caveat: DeepVariant models are trained and evaluated heavily on GIAB reference samples.
Performance on underrepresented populations or complex structural variant regions may
not match published benchmarks. Independent benchmarking on population-matched samples
is recommended before clinical deployment.

## Resource Requirements

| Data Type | RAM | CPU Time | GPU Time | Notes |
|-----------|-----|----------|----------|-------|
| WGS 30x | 64 GB | ~4-6 hrs | ~1-2 hrs | `--num_shards` scales make_examples linearly |
| WES | 32 GB | ~30 min | ~10 min | Much faster due to smaller target region |
| PacBio HiFi 30x | 64 GB | ~3-5 hrs | ~1-2 hrs | Fewer but longer reads |
| ONT 50x | 64 GB | ~6-8 hrs | ~2-3 hrs | Higher error rate means more candidate sites |

GPU acceleration primarily benefits the call_variants step. For large cohorts,
parallelizing across samples on CPU nodes may be more cost-effective than queuing
for GPU access.

## Complete Workflow Script

```bash
#!/bin/bash
set -euo pipefail

BAM=$1
REFERENCE=$2
OUTPUT_PREFIX=$3
MODEL_TYPE=${4:-WGS}
THREADS=${5:-16}

echo "=== DeepVariant: ${MODEL_TYPE} mode ==="

docker run -v "${PWD}:/data" google/deepvariant:1.6.1 \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=${MODEL_TYPE} \
    --ref=/data/${REFERENCE} \
    --reads=/data/${BAM} \
    --output_vcf=/data/${OUTPUT_PREFIX}.vcf.gz \
    --output_gvcf=/data/${OUTPUT_PREFIX}.g.vcf.gz \
    --intermediate_results_dir=/data/${OUTPUT_PREFIX}_tmp \
    --num_shards=${THREADS}

echo "=== Indexing ==="
bcftools index -t ${OUTPUT_PREFIX}.vcf.gz
bcftools index -t ${OUTPUT_PREFIX}.g.vcf.gz

echo "=== Statistics ==="
bcftools stats ${OUTPUT_PREFIX}.vcf.gz > ${OUTPUT_PREFIX}_stats.txt

echo "=== Complete ==="
echo "VCF: ${OUTPUT_PREFIX}.vcf.gz"
echo "gVCF: ${OUTPUT_PREFIX}.g.vcf.gz"
```

## Related Skills

- variant-calling/gatk-variant-calling - GATK alternative with joint calling ecosystem and VQSR integration
- variant-calling/variant-calling - bcftools calling for quick, lightweight analysis
- variant-calling/filtering-best-practices - Post-calling filtering strategies
- variant-calling/joint-calling - GATK joint genotyping alternative to GLnexus
- long-read-sequencing/clair3-variants - Long-read variant calling alternative, especially for ONT
