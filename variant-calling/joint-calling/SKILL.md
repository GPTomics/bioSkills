---
name: bio-variant-calling-joint-calling
description: Joint genotype calling across multiple samples using GATK CombineGVCFs and GenotypeGVCFs. Essential for cohort studies, population genetics, and leveraging VQSR. Use when performing joint genotyping across multiple samples.
tool_type: cli
primary_tool: GATK
---

## Version Compatibility

Reference examples tested with: GATK 4.5+, bcftools 1.19+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Joint Calling

**"Joint genotype my cohort samples"** → Combine per-sample gVCFs into a single cohort callset with consistent genotyping across all sites, enabling VQSR and population-level analysis.
- CLI: `gatk HaplotypeCaller -ERC GVCF` → `gatk GenomicsDBImport` → `gatk GenotypeGVCFs`

## Why Joint Calling Matters

Single-sample calling discards cross-sample evidence that is critical for accurate genotyping:

- **Statistical power from shared evidence** - A site with 2 alt reads in one sample is borderline and would typically be missed. If 50 other samples in the cohort also show 2 alt reads at that site, the evidence is overwhelming and the variant is clearly real. Joint calling aggregates this weak-per-sample signal into strong cohort-level evidence.
- **Genotype refinement via cohort priors** - Individual genotype likelihoods are combined with cohort allele frequencies as a Bayesian prior. A heterozygous call at a common variant site (AF=0.3) receives more support than the same call at a site with no other carriers. This prior dramatically improves accuracy for low-coverage samples.
- **Consistent site representation** - All samples are genotyped at the same sites, producing homozygous-reference calls where applicable. Without joint calling, a missing genotype is ambiguous: it could mean homozygous-reference or simply insufficient coverage. This "missing = reference" assumption is a common source of false negatives in downstream analysis.
- **VQSR eligibility** - Variant quality score recalibration requires cohort-level variant distributions and generally needs 30+ samples to build reliable models.

## Cohort Size Decision Table

| Cohort Size | Approach | Notes |
|---|---|---|
| <100 | CombineGVCFs or GenomicsDB | Either works; CombineGVCFs is simpler to manage |
| 100-10,000 | GenomicsDB + GenotypeGVCFs | Standard GATK Best Practices; shard by chromosome |
| 10,000-100,000 | GATK Biggest Practices | Heavily sharded and parallelized across intervals |
| >100,000 | DeepVariant + GLnexus, or Hail VDS | GATK becomes unwieldy at this scale; purpose-built tools required |

## Workflow Overview

```
Sample BAMs
    │
    ├── HaplotypeCaller (per-sample, -ERC GVCF)
    │   └── sample1.g.vcf.gz, sample2.g.vcf.gz, ...
    │
    ├── CombineGVCFs or GenomicsDBImport
    │   └── Combine into cohort database
    │
    ├── GenotypeGVCFs
    │   └── Joint genotyping
    │
    └── VQSR or Hard Filtering
        └── Final VCF
```

## Step 1: Per-Sample gVCF Generation

```bash
# Generate gVCF for each sample
gatk HaplotypeCaller \
    -R reference.fa \
    -I sample1.bam \
    -O sample1.g.vcf.gz \
    -ERC GVCF

# With intervals (faster)
gatk HaplotypeCaller \
    -R reference.fa \
    -I sample1.bam \
    -O sample1.g.vcf.gz \
    -ERC GVCF \
    -L intervals.bed
```

### Batch Processing

```bash
# Process all samples
for bam in *.bam; do
    sample=$(basename $bam .bam)
    gatk HaplotypeCaller \
        -R reference.fa \
        -I $bam \
        -O ${sample}.g.vcf.gz \
        -ERC GVCF &
done
wait
```

## Step 2a: CombineGVCFs (Small Cohorts)

For <100 samples:

```bash
gatk CombineGVCFs \
    -R reference.fa \
    -V sample1.g.vcf.gz \
    -V sample2.g.vcf.gz \
    -V sample3.g.vcf.gz \
    -O cohort.g.vcf.gz
```

### From Sample Map

```bash
# Create sample map file
# sample1    /path/to/sample1.g.vcf.gz
# sample2    /path/to/sample2.g.vcf.gz

ls *.g.vcf.gz | while read f; do
    echo -e "$(basename $f .g.vcf.gz)\t$f"
done > sample_map.txt

# Combine with -V for each
gatk CombineGVCFs \
    -R reference.fa \
    $(cat sample_map.txt | cut -f2 | sed 's/^/-V /') \
    -O cohort.g.vcf.gz
```

## Step 2b: GenomicsDBImport (Large Cohorts)

For >100 samples, use GenomicsDB:

```bash
# Create sample map
ls *.g.vcf.gz | while read f; do
    echo -e "$(basename $f .g.vcf.gz)\t$f"
done > sample_map.txt

# Import to GenomicsDB (per chromosome for parallelism)
gatk GenomicsDBImport \
    --sample-name-map sample_map.txt \
    --genomicsdb-workspace-path genomicsdb_chr1 \
    -L chr1 \
    --reader-threads 4

# Or all chromosomes
for chr in {1..22} X Y; do
    gatk GenomicsDBImport \
        --sample-name-map sample_map.txt \
        --genomicsdb-workspace-path genomicsdb_chr${chr} \
        -L chr${chr} &
done
wait
```

### Update GenomicsDB with New Samples

```bash
gatk GenomicsDBImport \
    --genomicsdb-update-workspace-path genomicsdb_chr1 \
    --sample-name-map new_samples.txt \
    -L chr1
```

### GenomicsDB Critical Caveats

GenomicsDB is powerful but has sharp edges that can cause data loss or silent failures:

- **No sample replacement** - Existing samples cannot be updated or overwritten. Only new samples with different names can be added. To fix a sample, the entire workspace must be recreated.
- **Intervals locked at import time** - The genomic intervals specified during the initial import cannot be changed on incremental updates. Adding new regions requires reimporting from scratch.
- **Fragment accumulation** - Each incremental batch creates a new database fragment. After thousands of incremental additions, file handle exhaustion becomes likely. Run `--consolidate` periodically to merge fragments.
- **Corruption risk on failed adds** - A failed incremental import can leave the datastore in an inconsistent state. Always backup the workspace directory before running `--genomicsdb-update-workspace-path`.
- **Batch size for memory** - Use `--batch-size 50` (the default) to control memory consumption. Larger batches load more gVCFs simultaneously and can exhaust heap space.

## Step 3: GenotypeGVCFs

### From Combined gVCF

```bash
gatk GenotypeGVCFs \
    -R reference.fa \
    -V cohort.g.vcf.gz \
    -O cohort.vcf.gz
```

### From GenomicsDB

```bash
gatk GenotypeGVCFs \
    -R reference.fa \
    -V gendb://genomicsdb_chr1 \
    -O chr1.vcf.gz

# All chromosomes
for chr in {1..22} X Y; do
    gatk GenotypeGVCFs \
        -R reference.fa \
        -V gendb://genomicsdb_chr${chr} \
        -O chr${chr}.vcf.gz &
done
wait

# Merge chromosomes
bcftools concat chr{1..22}.vcf.gz chrX.vcf.gz chrY.vcf.gz \
    -Oz -o cohort.vcf.gz
```

### With Allele-Specific Annotations

For larger cohorts where multiallelic sites are common, allele-specific annotations allow VQSR to evaluate each allele independently rather than penalizing a good allele because a co-occurring allele is poor:

```bash
gatk GenotypeGVCFs \
    -R reference.fa \
    -V gendb://genomicsdb \
    -O cohort.vcf.gz \
    -G StandardAnnotation \
    -G AS_StandardAnnotation
```

When allele-specific annotations are present, use `-AS` mode in VariantRecalibrator and ApplyVQSR for allele-level filtering.

## Step 4: Filtering

### VQSR (Recommended for >30 Samples)

```bash
# SNPs
gatk VariantRecalibrator \
    -R reference.fa \
    -V cohort.vcf.gz \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap.vcf.gz \
    --resource:omni,known=false,training=true,truth=false,prior=12.0 omni.vcf.gz \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp.vcf.gz \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode SNP \
    -O snps.recal \
    --tranches-file snps.tranches

gatk ApplyVQSR \
    -R reference.fa \
    -V cohort.vcf.gz \
    --recal-file snps.recal \
    --tranches-file snps.tranches \
    -mode SNP \
    --truth-sensitivity-filter-level 99.5 \
    -O cohort.snps.vcf.gz

# Indels
gatk VariantRecalibrator \
    -R reference.fa \
    -V cohort.snps.vcf.gz \
    --resource:mills,known=false,training=true,truth=true,prior=12.0 mills.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp.vcf.gz \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode INDEL \
    -O indels.recal \
    --tranches-file indels.tranches

gatk ApplyVQSR \
    -R reference.fa \
    -V cohort.snps.vcf.gz \
    --recal-file indels.recal \
    --tranches-file indels.tranches \
    -mode INDEL \
    --truth-sensitivity-filter-level 99.0 \
    -O cohort.filtered.vcf.gz
```

### Hard Filtering (Small Cohorts)

```bash
# See filtering-best-practices skill
gatk VariantFiltration \
    -R reference.fa \
    -V cohort.vcf.gz \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "FS > 60.0" --filter-name "FS60" \
    --filter-expression "MQ < 40.0" --filter-name "MQ40" \
    -O cohort.filtered.vcf.gz
```

## Batch Effects in Joint Calling

Joint genotyping mitigates most batch effects because it re-evaluates genotype likelihoods across all samples simultaneously, recalibrating quality scores against the full cohort distribution. However, certain batch effects persist through joint calling because they affect the underlying read data, not the genotyping model:

- **Different library prep protocols** - PCR-free vs PCR-based libraries produce different duplicate and error profiles
- **Different capture kits (WES)** - Exome kits target different regions; sites outside the intersection have systematically missing data in some batches
- **Significantly different coverage distributions** - 10x WGS samples mixed with 30x samples will have systematically different genotype quality at heterozygous sites
- **Different reference genome versions** - Mixing GRCh37 and GRCh38 alignments is not valid; all samples must use the same reference
- **Mixing WGS and WES** - Fundamentally different coverage profiles; off-target WES regions behave like very-low-coverage WGS

Mitigation: process all samples through an identical upstream pipeline (same aligner, same duplicate marking, same BQSR resources). If batches are unavoidable, include batch as a covariate in downstream association or differential analyses.

## When to Re-genotype

| Scenario | Action | Rationale |
|---|---|---|
| Adding new samples | Re-genotype (GenomicsDB incremental add + GenotypeGVCFs on full database) | New samples change cohort allele frequencies, improving all genotype calls |
| Changing reference genome | Full reprocess from alignment | gVCF coordinates are reference-specific |
| Updating caller version | Optional but recommended for consistency | Different caller versions may produce different quality scores; mixing versions adds noise |
| Adding new genomic intervals | Reimport from scratch | GenomicsDB intervals are locked at initial import; incremental update cannot expand them |

## DeepVariant + GLnexus Alternative

For cohorts using DeepVariant as the primary caller, GLnexus provides a scalable joint calling pathway that has been validated at 250,000+ samples:

```bash
# Step 1: Run DeepVariant per sample to produce gVCFs
deepvariant --model_type=WGS \
    --ref=reference.fa --reads=sample.bam \
    --output_vcf=sample.vcf.gz --output_gvcf=sample.g.vcf.gz

# Step 2: Joint call with GLnexus (pre-tuned configs available)
glnexus_cli --config DeepVariantWGS --bed intervals.bed \
    sample1.g.vcf.gz sample2.g.vcf.gz ... | bcftools view - | bgzip -c > cohort.vcf.gz
```

Pre-tuned GLnexus configs: `DeepVariantWGS` for whole-genome, `DeepVariantWES` for whole-exome. These configs encode appropriate genotype quality thresholds and multiallelic handling tuned for DeepVariant output.

## Complete Pipeline Script

**Goal:** Run the full joint calling workflow from BAMs to filtered cohort VCF.

**Approach:** Generate per-sample gVCFs, import into GenomicsDB, joint genotype, then index and compute statistics.

```bash
#!/bin/bash
set -euo pipefail

REFERENCE=$1
OUTPUT_DIR=$2
THREADS=16

mkdir -p $OUTPUT_DIR/{gvcfs,genomicsdb,vcfs}

echo "=== Step 1: Generate gVCFs ==="
for bam in data/*.bam; do
    sample=$(basename $bam .bam)
    gatk HaplotypeCaller \
        -R $REFERENCE \
        -I $bam \
        -O $OUTPUT_DIR/gvcfs/${sample}.g.vcf.gz \
        -ERC GVCF &

    # Limit parallelism
    while [ $(jobs -r | wc -l) -ge $THREADS ]; do sleep 1; done
done
wait

echo "=== Step 2: Create sample map ==="
ls $OUTPUT_DIR/gvcfs/*.g.vcf.gz | while read f; do
    echo -e "$(basename $f .g.vcf.gz)\t$(realpath $f)"
done > $OUTPUT_DIR/sample_map.txt

echo "=== Step 3: GenomicsDBImport ==="
gatk GenomicsDBImport \
    --sample-name-map $OUTPUT_DIR/sample_map.txt \
    --genomicsdb-workspace-path $OUTPUT_DIR/genomicsdb \
    -L intervals.bed \
    --reader-threads 4

echo "=== Step 4: Joint genotyping ==="
gatk GenotypeGVCFs \
    -R $REFERENCE \
    -V gendb://$OUTPUT_DIR/genomicsdb \
    -O $OUTPUT_DIR/vcfs/cohort.vcf.gz

echo "=== Step 5: Index ==="
bcftools index -t $OUTPUT_DIR/vcfs/cohort.vcf.gz

echo "=== Statistics ==="
bcftools stats $OUTPUT_DIR/vcfs/cohort.vcf.gz > $OUTPUT_DIR/vcfs/cohort_stats.txt

echo "=== Complete ==="
echo "Joint VCF: $OUTPUT_DIR/vcfs/cohort.vcf.gz"
```

## Tips

### Memory for Large Cohorts

```bash
# Increase Java heap for GenotypeGVCFs (default 4g is often insufficient for >500 samples)
gatk --java-options "-Xmx64g" GenotypeGVCFs ...

# For GenomicsDBImport, --batch-size controls how many gVCFs are loaded simultaneously
gatk GenomicsDBImport --batch-size 50 ...
```

## Related Skills

- variant-calling/gatk-variant-calling - Single-sample HaplotypeCaller workflow
- variant-calling/deepvariant - DeepVariant caller for use with GLnexus pathway
- variant-calling/filtering-best-practices - VQSR and hard filtering strategies
- variant-calling/vcf-manipulation - Merging and subsetting joint-called VCFs
- population-genetics/plink-basics - Population analysis of joint calls
- workflows/fastq-to-variants - End-to-end germline pipeline
