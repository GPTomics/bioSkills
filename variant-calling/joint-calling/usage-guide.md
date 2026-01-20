# Joint Calling Usage Guide

Multi-sample joint genotyping for improved variant detection.

## Prerequisites

```bash
conda install -c bioconda gatk4
```

## Why Joint Calling?

- **Better sensitivity**: Leverage information across samples
- **Consistent sites**: Same positions called in all samples
- **VQSR eligible**: Machine learning filtering requires cohorts
- **Population frequencies**: Calculate allele frequencies across cohort

## Workflow Overview

```
Sample BAMs → HaplotypeCaller (GVCF mode) → CombineGVCFs → GenotypeGVCFs → Cohort VCF
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
| < 50 | CombineGVCFs | Simple, single command |
| 50-1000 | GenomicsDBImport | Scalable, interval-based |
| > 1000 | GenomicsDBImport + batches | Process in batches |

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

## See Also

- [GATK joint calling tutorial](https://gatk.broadinstitute.org/hc/en-us/articles/360035531152)
