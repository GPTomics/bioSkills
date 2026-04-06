---
name: bio-gatk-variant-calling
description: Variant calling with GATK HaplotypeCaller following best practices. Covers germline SNP/indel calling, GVCF workflow for cohorts, joint genotyping, and variant quality score recalibration (VQSR). Use when calling variants with GATK HaplotypeCaller.
tool_type: cli
primary_tool: gatk
---

## Version Compatibility

Reference examples tested with: GATK 4.5+, bcftools 1.19+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# GATK Variant Calling

GATK HaplotypeCaller performs local de novo assembly of haplotypes in active regions, making it more accurate than pileup-based callers (bcftools) for indels and complex variants.

## Pipeline Decision Tree

```
What is the analysis context?
├── Single sample, highest accuracy → DRAGEN-GATK mode (hard filter on QUAL)
├── Cohort <1000, human → Standard Best Practices (GVCF → joint genotype → VQSR)
├── Cohort >1000, human → GATK "Biggest Practices" or DeepVariant + GLnexus
├── Non-human organism → Hard filtering (no VQSR training resources)
├── Targeted panel / small exome → Hard filtering (too few variants for VQSR)
└── Somatic variants → Mutect2 (not HaplotypeCaller)
```

## DRAGEN-GATK Mode (Current Recommendation)

DRAGEN-GATK mode is the recommended approach for single-sample germline calling. It incorporates three innovations that improve accuracy:
- **BQD (Base Quality Dropout)**: Corrects systematic base quality errors at read ends
- **FRD (Foreign Read Detection)**: Identifies contaminating reads at the genotyping level
- **DragSTR**: Improved short tandem repeat handling with a Markov model

Key differences from standard mode: BQSR is no longer needed (error correction happens during calling), and QUAL scores are more calibrated, making simple hard filtering on QUAL sufficient without VQSR.

```bash
gatk HaplotypeCaller \
    -R reference.fa \
    -I sample.bam \
    -O sample.vcf.gz \
    --dragen-mode
```

For DRAGEN-GVCF mode (cohort workflows):
```bash
gatk HaplotypeCaller \
    -R reference.fa \
    -I sample.bam \
    -O sample.g.vcf.gz \
    -ERC GVCF \
    --dragen-mode
```

Note: DRAGEN mode is available only for single-sample calling. Joint genotyping with GenotypeGVCFs does not use DRAGEN mode but benefits from the improved per-sample GVCFs.

## Prerequisites

BAM files should be preprocessed:
1. Mark duplicates (always required)
2. Base quality score recalibration (BQSR) -- recommended for standard mode, NOT needed for DRAGEN-GATK mode

## Single-Sample Calling

**Goal:** Call germline SNPs and indels from a single sample using HaplotypeCaller.

**Approach:** Run local de novo assembly of haplotypes in active regions to detect variants with optional annotation enrichment.

**"Call variants from my BAM file using GATK"** → Perform local haplotype assembly and genotyping on aligned reads using HaplotypeCaller.

### Basic HaplotypeCaller

```bash
gatk HaplotypeCaller \
    -R reference.fa \
    -I sample.bam \
    -O sample.vcf.gz
```

### With Standard Annotations

```bash
gatk HaplotypeCaller \
    -R reference.fa \
    -I sample.bam \
    -O sample.vcf.gz \
    -A Coverage \
    -A QualByDepth \
    -A FisherStrand \
    -A StrandOddsRatio \
    -A MappingQualityRankSumTest \
    -A ReadPosRankSumTest
```

### Target Intervals (Exome/Panel)

```bash
gatk HaplotypeCaller \
    -R reference.fa \
    -I sample.bam \
    -L targets.interval_list \
    -O sample.vcf.gz
```

### Adjust Calling Confidence

```bash
gatk HaplotypeCaller \
    -R reference.fa \
    -I sample.bam \
    -O sample.vcf.gz \
    --standard-min-confidence-threshold-for-calling 20
```

## GVCF Workflow (Recommended for Cohorts)

**Goal:** Enable joint genotyping across a cohort by generating per-sample genomic VCFs.

**Approach:** Call each sample in GVCF mode (-ERC GVCF), combine into a GenomicsDB or merged GVCF, then jointly genotype.

The GVCF workflow enables joint genotyping across samples for better variant calls.

### Step 1: Generate GVCFs per Sample

```bash
gatk HaplotypeCaller \
    -R reference.fa \
    -I sample.bam \
    -O sample.g.vcf.gz \
    -ERC GVCF
```

### Step 2: Combine GVCFs (GenomicsDBImport)

```bash
# Create sample map file
# sample_map.txt:
# sample1    /path/to/sample1.g.vcf.gz
# sample2    /path/to/sample2.g.vcf.gz

gatk GenomicsDBImport \
    --genomicsdb-workspace-path genomicsdb \
    --sample-name-map sample_map.txt \
    -L intervals.interval_list
```

### Alternative: CombineGVCFs (smaller cohorts)

```bash
gatk CombineGVCFs \
    -R reference.fa \
    -V sample1.g.vcf.gz \
    -V sample2.g.vcf.gz \
    -V sample3.g.vcf.gz \
    -O cohort.g.vcf.gz
```

### Step 3: Joint Genotyping

```bash
# From GenomicsDB
gatk GenotypeGVCFs \
    -R reference.fa \
    -V gendb://genomicsdb \
    -O cohort.vcf.gz

# From combined GVCF
gatk GenotypeGVCFs \
    -R reference.fa \
    -V cohort.g.vcf.gz \
    -O cohort.vcf.gz
```

## Variant Quality Score Recalibration (VQSR)

**Goal:** Apply machine learning-based variant filtering using known truth/training sets.

**Approach:** Build a Gaussian mixture model from annotation values at known sites, then apply a sensitivity threshold to classify variants.

### When to Use VQSR vs Hard Filtering

| Context | Filtering Method | Rationale |
|---------|-----------------|-----------|
| Human WGS cohort, >30 samples | VQSR | Enough variants for model training; truth sets available |
| Human WGS, single sample, DRAGEN mode | Hard filter on QUAL | DRAGEN QUAL scores are well-calibrated |
| Human exome, >30 samples | VQSR (with `--max-gaussians 4` for indels) | Fewer variants but usually sufficient |
| <30 exomes or gene panels | Hard filtering | Too few variants for reliable GMM training |
| Non-model organisms | Hard filtering | No HapMap/1000G training resources |
| Somatic calling (Mutect2) | FilterMutectCalls | Dedicated somatic filtering, not VQSR |

Sensitivity tranche selection: 99.5% for SNPs (0.5% true positives may be filtered) and 99.0% for indels (indels have higher error rates). Adjust based on the analysis goal -- use 99.7% for discovery, 99.0% for clinical stringency.

### SNP Recalibration

```bash
# Build SNP model
# Annotations: QD, MQ are most informative; FS/SOR detect strand bias; RankSum tests detect systematic ref/alt differences
gatk VariantRecalibrator \
    -R reference.fa \
    -V cohort.vcf.gz \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap.vcf.gz \
    --resource:omni,known=false,training=true,truth=false,prior=12.0 omni.vcf.gz \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp.vcf.gz \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode SNP \
    -O snp.recal \
    --tranches-file snp.tranches

# Apply SNP filter (99.5% sensitivity retains 99.5% of true positives)
gatk ApplyVQSR \
    -R reference.fa \
    -V cohort.vcf.gz \
    -O cohort.snp_recal.vcf.gz \
    --recal-file snp.recal \
    --tranches-file snp.tranches \
    --truth-sensitivity-filter-level 99.5 \
    -mode SNP
```

### Allele-Specific VQSR (Recommended for Large Cohorts)

Standard VQSR evaluates all alleles at a site together. Allele-specific (AS) filtering evaluates each allele independently -- critical at multiallelic sites where a true variant and an artifact may co-occur. Increasingly important as cohort size grows (more multiallelic sites).

```bash
# Requires AS annotations from GenotypeGVCFs (-G AS_StandardAnnotation)
gatk VariantRecalibrator \
    -R reference.fa \
    -V cohort.vcf.gz \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap.vcf.gz \
    --resource:omni,known=false,training=true,truth=false,prior=12.0 omni.vcf.gz \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp.vcf.gz \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode SNP -AS \
    -O snp.as.recal \
    --tranches-file snp.as.tranches

gatk ApplyVQSR \
    -R reference.fa \
    -V cohort.vcf.gz \
    -O cohort.snp_as_recal.vcf.gz \
    --recal-file snp.as.recal \
    --tranches-file snp.as.tranches \
    --truth-sensitivity-filter-level 99.5 \
    -mode SNP -AS
```

### Indel Recalibration

```bash
# Build Indel model
gatk VariantRecalibrator \
    -R reference.fa \
    -V cohort.snp_recal.vcf.gz \
    --resource:mills,known=false,training=true,truth=true,prior=12.0 Mills.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp.vcf.gz \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode INDEL \
    --max-gaussians 4 \
    -O indel.recal \
    --tranches-file indel.tranches

# Apply Indel filter
gatk ApplyVQSR \
    -R reference.fa \
    -V cohort.snp_recal.vcf.gz \
    -O cohort.vqsr.vcf.gz \
    --recal-file indel.recal \
    --tranches-file indel.tranches \
    --truth-sensitivity-filter-level 99.0 \
    -mode INDEL
```

## Hard Filtering (When VQSR Not Suitable)

**Goal:** Apply fixed-threshold filters when the dataset is too small for VQSR.

**Approach:** Separate SNPs and indels (they have fundamentally different annotation distributions), apply GATK-recommended annotation thresholds, then merge results.

For small datasets, exomes, single samples, non-model organisms, or DRAGEN-mode output. These thresholds are starting points -- always examine annotation distributions in the specific dataset and adjust accordingly.

### Extract SNPs and Indels

```bash
gatk SelectVariants \
    -R reference.fa \
    -V cohort.vcf.gz \
    --select-type-to-include SNP \
    -O snps.vcf.gz

gatk SelectVariants \
    -R reference.fa \
    -V cohort.vcf.gz \
    --select-type-to-include INDEL \
    -O indels.vcf.gz
```

### Apply Hard Filters

```bash
# Filter SNPs
gatk VariantFiltration \
    -R reference.fa \
    -V snps.vcf.gz \
    -O snps.filtered.vcf.gz \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "FS > 60.0" --filter-name "FS60" \
    --filter-expression "MQ < 40.0" --filter-name "MQ40" \
    --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    --filter-expression "SOR > 3.0" --filter-name "SOR3"

# Filter Indels
gatk VariantFiltration \
    -R reference.fa \
    -V indels.vcf.gz \
    -O indels.filtered.vcf.gz \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "FS > 200.0" --filter-name "FS200" \
    --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    --filter-expression "SOR > 10.0" --filter-name "SOR10"
```

### Merge Filtered Variants

```bash
gatk MergeVcfs \
    -I snps.filtered.vcf.gz \
    -I indels.filtered.vcf.gz \
    -O cohort.filtered.vcf.gz
```

## Base Quality Score Recalibration (BQSR)

Optional preprocessing for standard mode (NOT needed for DRAGEN-GATK mode). Corrects systematic base quality errors using known variant sites.

```bash
gatk BaseRecalibrator -R reference.fa -I sample.bam \
    --known-sites dbsnp.vcf.gz --known-sites known_indels.vcf.gz -O recal_data.table

gatk ApplyBQSR -R reference.fa -I sample.bam \
    --bqsr-recal-file recal_data.table -O sample.recal.bam
```

## Parallel Processing

### Scatter by Interval

```bash
for interval in chr{1..22} chrX chrY; do
    gatk HaplotypeCaller -R reference.fa -I sample.bam -L $interval \
        -O sample.${interval}.g.vcf.gz -ERC GVCF &
done
wait
gatk GatherVcfs $(for c in chr{1..22} chrX chrY; do echo "-I sample.${c}.g.vcf.gz"; done) -O sample.g.vcf.gz
```

### Native PairHMM Threads

```bash
gatk HaplotypeCaller -R reference.fa -I sample.bam -O sample.vcf.gz --native-pair-hmm-threads 4
```

## CNN Score Variant Filter (Deep Learning)

Alternative to VQSR using convolutional neural networks. Less established than VQSR but does not require cohort data or truth sets. NVScoreVariants (NVIDIA) is the successor with improved handling of depth variation.

```bash
gatk CNNScoreVariants -R reference.fa -V cohort.vcf.gz -O cohort.cnn_scored.vcf.gz --tensor-type reference

gatk FilterVariantTranches -V cohort.cnn_scored.vcf.gz -O cohort.cnn_filtered.vcf.gz \
    --resource hapmap.vcf.gz --resource mills.vcf.gz \
    --info-key CNN_1D --snp-tranche 99.95 --indel-tranche 99.4
```

## Complete Single-Sample Pipeline (DRAGEN Mode)

**Goal:** Run the current recommended GATK workflow for single-sample germline calling.

**Approach:** Mark duplicates, call with HaplotypeCaller in DRAGEN mode (no BQSR needed), genotype, and hard filter on QUAL.

```bash
#!/bin/bash
SAMPLE=$1
REF=reference.fa

# Mark duplicates (required for all modes)
gatk MarkDuplicates -I ${SAMPLE}.bam -O ${SAMPLE}.markdup.bam -M ${SAMPLE}.dup_metrics.txt
samtools index ${SAMPLE}.markdup.bam

# Call variants in DRAGEN mode (BQSR not needed)
gatk HaplotypeCaller -R $REF -I ${SAMPLE}.markdup.bam \
    -O ${SAMPLE}.g.vcf.gz -ERC GVCF --dragen-mode

# Single-sample genotyping
gatk GenotypeGVCFs -R $REF -V ${SAMPLE}.g.vcf.gz \
    -O ${SAMPLE}.vcf.gz

# DRAGEN mode: hard filter on QUAL is sufficient
gatk VariantFiltration -R $REF -V ${SAMPLE}.vcf.gz \
    -O ${SAMPLE}.filtered.vcf.gz \
    --filter-expression "QUAL < 10.4139" --filter-name "LowQual"
```

## Complete Single-Sample Pipeline (Standard Mode)

**Goal:** Run the classic GATK best practices workflow with BQSR.

**Approach:** Chain BaseRecalibrator, ApplyBQSR, HaplotypeCaller, GenotypeGVCFs, and hard filtering.

```bash
#!/bin/bash
SAMPLE=$1
REF=reference.fa
DBSNP=dbsnp.vcf.gz
KNOWN_INDELS=known_indels.vcf.gz

# BQSR
gatk BaseRecalibrator -R $REF -I ${SAMPLE}.bam \
    --known-sites $DBSNP --known-sites $KNOWN_INDELS \
    -O ${SAMPLE}.recal.table

gatk ApplyBQSR -R $REF -I ${SAMPLE}.bam \
    --bqsr-recal-file ${SAMPLE}.recal.table \
    -O ${SAMPLE}.recal.bam

# Call variants
gatk HaplotypeCaller -R $REF -I ${SAMPLE}.recal.bam \
    -O ${SAMPLE}.g.vcf.gz -ERC GVCF

# Single-sample genotyping
gatk GenotypeGVCFs -R $REF -V ${SAMPLE}.g.vcf.gz \
    -O ${SAMPLE}.vcf.gz

# Hard filter (SNPs and indels need different thresholds -- see filtering-best-practices)
gatk SelectVariants -R $REF -V ${SAMPLE}.vcf.gz --select-type-to-include SNP -O ${SAMPLE}.snps.vcf.gz
gatk SelectVariants -R $REF -V ${SAMPLE}.vcf.gz --select-type-to-include INDEL -O ${SAMPLE}.indels.vcf.gz

gatk VariantFiltration -R $REF -V ${SAMPLE}.snps.vcf.gz -O ${SAMPLE}.snps.filtered.vcf.gz \
    --filter-expression "QD < 2.0" --filter-name "LowQD" \
    --filter-expression "FS > 60.0" --filter-name "HighFS" \
    --filter-expression "MQ < 40.0" --filter-name "LowMQ" \
    --filter-expression "SOR > 3.0" --filter-name "HighSOR" \
    --filter-expression "MQRankSum < -12.5" --filter-name "LowMQRS" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "LowRPRS"

gatk VariantFiltration -R $REF -V ${SAMPLE}.indels.vcf.gz -O ${SAMPLE}.indels.filtered.vcf.gz \
    --filter-expression "QD < 2.0" --filter-name "LowQD" \
    --filter-expression "FS > 200.0" --filter-name "HighFS" \
    --filter-expression "SOR > 10.0" --filter-name "HighSOR" \
    --filter-expression "ReadPosRankSum < -20.0" --filter-name "LowRPRS"

gatk MergeVcfs -I ${SAMPLE}.snps.filtered.vcf.gz -I ${SAMPLE}.indels.filtered.vcf.gz \
    -O ${SAMPLE}.filtered.vcf.gz
```

## Key Annotations and Filter Rationale

| Annotation | SNP Threshold | Indel Threshold | What It Detects |
|-----------|---------------|-----------------|-----------------|
| QD | < 2.0 | < 2.0 | Low variant quality relative to depth (not supported by reads) |
| FS | > 60.0 | > 200.0 | Strand bias: variant reads predominantly on one strand (artifact). Indels naturally show more strand bias, hence higher threshold |
| SOR | > 3.0 | > 10.0 | Same as FS but handles high-depth sites better using symmetric odds ratio |
| MQ | < 40.0 | N/A | Reads map ambiguously (paralogous regions, segmental duplications) |
| MQRankSum | < -12.5 | < -12.5 | Alt-supporting reads have worse mapping quality than ref reads |
| ReadPosRankSum | < -8.0 | < -20.0 | Variant appears only at read ends (misalignment artifact). More permissive for indels because indel reads often cluster at specific positions |

These thresholds are conservative starting points. For specific datasets, plot annotation histograms for known true/false variants and adjust thresholds. Over-filtering disproportionately affects rare variants -- consider tiered filtering (stringent for common variants, permissive for rare variants) in clinical or discovery contexts.

## Resource Files

| Resource | Use |
|----------|-----|
| dbSNP | Known variants (prior=2.0) |
| HapMap | Training/truth SNPs (prior=15.0) |
| Omni | Training SNPs (prior=12.0) |
| 1000G SNPs | Training SNPs (prior=10.0) |
| Mills Indels | Training/truth indels (prior=12.0) |

## Related Skills

- variant-calling/variant-calling - bcftools alternative (faster, less accurate)
- variant-calling/deepvariant - Deep learning alternative (highest accuracy)
- variant-calling/joint-calling - Multi-sample joint genotyping workflow
- variant-calling/filtering-best-practices - Post-calling filtering details
- variant-calling/variant-normalization - Normalize before annotation
- variant-calling/variant-annotation - Annotate final calls with VEP/SnpEff
