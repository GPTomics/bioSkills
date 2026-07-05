# Joint Calling - Usage Guide

## Overview

Multi-sample joint genotyping combines information across samples for improved variant detection, consistent site calling, and eligibility for cohort-level filtering (VQSR/VETS). It is not the same operation as merging single-sample VCFs: joint genotyping borrows evidence across samples to rescue low-coverage het sites and produces a squared-off sample-by-site matrix where an unassessed cell (`./.`) is distinguishable from a confident hom-ref (`0/0`). The GVCF workflow also solves the N+1 problem, letting a new sample be added without re-calling the whole cohort.

## Prerequisites

```bash
conda install -c bioconda gatk4
```

## Quick Start

Tell your AI agent what you want to do:
- "Set up joint calling for my cohort of 50 samples"
- "Combine GVCFs and run GenotypeGVCFs"
- "Scale joint calling for 1000+ samples by chromosome"
- "Choose between CombineGVCFs and GenomicsDBImport"
- "Add a new sample to my cohort without re-calling everyone"
- "Joint call DeepVariant gVCFs with GLnexus instead of GATK"
- "Explain why I should not just bcftools merge my single-sample VCFs"

## Example Prompts

### Small Cohort
> "Combine these 10 sample GVCFs and run joint genotyping"

> "Run CombineGVCFs followed by GenotypeGVCFs for my family trio"

### Large Cohort
> "Set up GenomicsDBImport for 500 samples, processing by chromosome"

> "Create a sample map file and import into GenomicsDB"

### GVCF Generation
> "Generate GVCFs from all BAM files for joint calling"

> "Run HaplotypeCaller in GVCF mode for each sample in the cohort"

### Post-Processing
> "Apply VQSR to my joint-called cohort VCF"

> "Joint call my samples and then apply hard filters"

### Scaling and Alternatives
> "Shard GenomicsDBImport and GenotypeGVCFs by chromosome for my 5000-sample cohort"

> "Reblock my gVCFs and use GnarlyGenotyper for a biobank-scale cohort"

> "Run DeepVariant per sample and joint call with GLnexus DeepVariantWGS config"

## What the Agent Will Do

1. Generate per-sample GVCFs using HaplotypeCaller with -ERC GVCF
2. Select appropriate combining method (CombineGVCFs vs GenomicsDBImport)
3. Combine GVCFs across samples
4. Run GenotypeGVCFs to produce the cohort VCF
5. Apply VQSR or hard filtering based on cohort size
6. Validate output with variant statistics

## Tips

- Use GenomicsDBImport for cohorts >100 samples (locus-centric TileDB storage scales); CombineGVCFs is simpler for small cohorts
- Never `bcftools merge` single-sample VCFs into a cohort matrix: absent records become `./.` (missing), not `0/0` (hom-ref), corrupting allele frequencies
- Cap the JVM heap (`--java-options -Xmx`) at ~80-90% of RAM for GenomicsDBImport; its heavy work is native C/C++, and an over-large heap triggers a native out-of-memory error
- Process large cohorts by interval (per-chromosome) to parallelize and manage memory; use a `--sample-name-map` file rather than thousands of `-V` arguments
- Always generate GVCFs (not VCFs) when planning joint calling; the gVCF is the reusable intermediate that solves the N+1 problem
- Expect `*` spanning-deletion alleles and recomputed GQ/PL in the joint VCF; special-case `*` before annotation
- VQSR/VETS need cohort-scale variants to train -- a single deep WGS supplies enough, but variant-poor exomes/panels need ~30+ samples; otherwise hard-filter
- For DeepVariant cohorts, GLnexus merges ~8x faster with ~7x smaller gVCFs than routing through GenotypeGVCFs (Yun 2020)

## Why Joint Calling?

- Better sensitivity: leverage information across samples
- Consistent sites: same positions called in all samples
- VQSR eligible: machine learning filtering requires cohorts
- Population frequencies: calculate allele frequencies across cohort

## Workflow Overview

```
Sample BAMs -> HaplotypeCaller (GVCF mode) -> CombineGVCFs -> GenotypeGVCFs -> Cohort VCF
```

## Step-by-Step

### 1. Generate per-sample GVCFs

```bash
gatk HaplotypeCaller \
    -R reference.fa \
    -I sample1.bam \
    -O sample1.g.vcf.gz \
    -ERC GVCF
```

### 2. Combine GVCFs

#### Option A: CombineGVCFs (small cohorts)

```bash
gatk CombineGVCFs \
    -R reference.fa \
    -V sample1.g.vcf.gz \
    -V sample2.g.vcf.gz \
    -V sample3.g.vcf.gz \
    -O combined.g.vcf.gz
```

#### Option B: GenomicsDBImport (large cohorts)

```bash
gatk GenomicsDBImport \
    -V sample1.g.vcf.gz \
    -V sample2.g.vcf.gz \
    -V sample3.g.vcf.gz \
    --genomicsdb-workspace-path genomicsdb \
    -L intervals.bed
```

### 3. Joint Genotyping

```bash
# From CombineGVCFs
gatk GenotypeGVCFs \
    -R reference.fa \
    -V combined.g.vcf.gz \
    -O cohort.vcf.gz

# From GenomicsDB
gatk GenotypeGVCFs \
    -R reference.fa \
    -V gendb://genomicsdb \
    -O cohort.vcf.gz
```

## Scaling Tips

| Cohort Size | Method | Notes |
|-------------|--------|-------|
| < 100 | CombineGVCFs | Simple, single command; portable plain gVCF |
| 100-10,000 | GenomicsDBImport | Scalable, interval-based; shard by chromosome |
| > 10,000 | ReblockGVCF + GnarlyGenotyper, or DeepVariant + GLnexus | Naive GenotypeGVCFs stops scaling |

### Large Cohort Strategy

```bash
# Import by chromosome
for chr in {1..22} X Y; do
    gatk GenomicsDBImport \
        --sample-name-map samples.map \
        --genomicsdb-workspace-path genomicsdb_chr${chr} \
        -L chr${chr}
done

# Genotype by chromosome, then merge
```

## Post-Processing

```bash
# Apply VQSR (requires large cohorts)
gatk VariantRecalibrator ...
gatk ApplyVQSR ...

# Or hard filter for small cohorts
gatk VariantFiltration \
    -V cohort.vcf.gz \
    --filter-expression "QD < 2.0" --filter-name "LowQD" \
    -O filtered.vcf.gz
```

## Related Skills

- variant-calling/gatk-variant-calling - Single-sample calling and per-sample gVCF generation
- variant-calling/filtering-best-practices - VQSR/VETS and hard filtering details
- variant-calling/deepvariant - DeepVariant + GLnexus alternative joint calling
- variant-calling/vcf-manipulation - Merge/subset semantics vs joint genotyping
- variant-calling/vcf-basics - Genotype-field grammar (`./.` vs `0/0`)
- population-genetics/plink-basics - Population analysis of joint calls
