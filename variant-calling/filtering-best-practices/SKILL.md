---
name: bio-variant-calling-filtering-best-practices
description: Comprehensive variant filtering including GATK VQSR, hard filters, bcftools expressions, and quality metric interpretation for SNPs and indels. Use when filtering variants using GATK best practices.
tool_type: mixed
primary_tool: bcftools
---

## Version Compatibility

Reference examples tested with: GATK 4.5+, bcftools 1.19+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Variant Filtering Best Practices

## Filter Selection Decision Tree

```
Is this somatic data?
├── Yes → FilterMutectCalls (GATK), not VQSR or hard filters
└── No (germline) →
    Was DRAGEN-GATK mode used?
    ├── Yes → Hard filter on QUAL only (DRAGEN QUAL is well-calibrated)
    └── No →
        Is dataset large enough for VQSR?
        ├── Yes (>30 WGS/exomes, human, truth sets available) → VQSR
        │   └── Large cohort? → Use allele-specific VQSR (-AS flag)
        └── No → Hard filtering
            ├── Non-model organism → Hard filtering only (no training resources)
            └── Targeted panel → Hard filtering (too few variants for VQSR model)
```

VQSR requires a Gaussian mixture model trained on known truth sets (HapMap, 1000G, dbSNP).
With fewer than ~30 samples, the model cannot learn the annotation distributions reliably.
Targeted panels produce too few variants for stable model convergence, regardless of sample count.

## GATK Hard Filter Thresholds

**Goal:** Apply GATK-recommended annotation thresholds to separate true variants from artifacts.

**Approach:** Use VariantFiltration with per-metric filter expressions for SNPs and indels separately.

**"Filter my variants using GATK best practices"** → Apply fixed annotation thresholds (QD, FS, MQ, SOR, RankSum) to flag low-quality variants.

```bash
# SNPs
gatk VariantFiltration \
    -R reference.fa \
    -V raw_snps.vcf \
    -O filtered_snps.vcf \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "FS > 60.0" --filter-name "FS60" \
    --filter-expression "MQ < 40.0" --filter-name "MQ40" \
    --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    --filter-expression "SOR > 3.0" --filter-name "SOR3"

# Indels
gatk VariantFiltration \
    -R reference.fa \
    -V raw_indels.vcf \
    -O filtered_indels.vcf \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "FS > 200.0" --filter-name "FS200" \
    --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    --filter-expression "SOR > 10.0" --filter-name "SOR10"
```

## Understanding Quality Metrics

| Metric | Threshold | Rationale |
|--------|-----------|-----------|
| QD | <2.0 | Quality normalized by depth; preferred over raw QUAL because high-depth sites inflate QUAL regardless of evidence quality. QD < 2.0 means variant quality is not supported by sufficient per-read evidence. |
| FS | >60 (SNP), >200 (indel) | Fisher strand bias p-value (Phred-scaled). Indels naturally exhibit more strand bias due to alignment artifacts near insertions/deletions, hence the relaxed indel threshold. |
| SOR | >3.0 (SNP), >10.0 (indel) | Symmetric odds ratio for strand bias. Handles high-depth sites better than FS by avoiding the inflated p-values that Fisher's exact test produces at high counts. |
| MQ | <40.0 | RMS mapping quality across all reads. Low MQ suggests reads map ambiguously, indicating paralogous regions or segmental duplications. |
| MQRankSum | <-12.5 | Compares mapping quality of alt-supporting vs ref-supporting reads. Large negative values mean alt reads map significantly worse, suggesting they originate from mismapped paralogs. |
| ReadPosRankSum | <-8.0 (SNP), <-20.0 (indel) | Position within read where the variant allele appears. Large negative values mean the variant clusters at read ends, a hallmark of misalignment or sequencing error. Indels tolerate a wider range because alignment uncertainty near indels shifts read positions. |
| DP | Sample-specific | Extreme depth (>2x mean or <0.3x mean) suggests collapsed repeats or poor capture. Filtering on depth alone removes real variants in duplicated regions; always combine with other annotations. |
| GQ | <20 | Phred-scaled confidence in the assigned genotype. GQ 20 = 99% confidence. |

## bcftools filter

**Goal:** Filter variants using bcftools expression syntax with soft or hard removal.

**Approach:** Use -e (exclude) or -i (include) with expressions on QUAL, INFO, and FORMAT fields; use -s for soft filtering.

### Soft vs Hard Filtering

```bash
# Hard filter (remove variants)
bcftools filter -e 'QUAL<30' input.vcf.gz -o filtered.vcf

# Soft filter (mark, don't remove)
bcftools filter -s 'LowQual' -e 'QUAL<30' input.vcf.gz -o marked.vcf
# Variants failing filter get "LowQual" in FILTER column

# Include instead of exclude
bcftools filter -i 'QUAL>=30' input.vcf.gz -o filtered.vcf
```

### Expression Syntax

| Operator | Meaning |
|----------|---------|
| `<`, `<=`, `>`, `>=` | Comparison |
| `=`, `==` | Equals |
| `!=` | Not equals |
| `&&`, `\|\|` | AND, OR |
| `!` | NOT |

### Aggregate Functions

| Function | Description |
|----------|-------------|
| `MIN(x)` | Minimum across samples |
| `MAX(x)` | Maximum across samples |
| `AVG(x)` | Average across samples |
| `SUM(x)` | Sum across samples |

### Common bcftools Filters

```bash
# Basic quality filter
bcftools filter -i 'QUAL>30 && DP>10' input.vcf -o filtered.vcf

# Complex filter with multiple metrics
bcftools filter -i 'QUAL>30 && INFO/DP>10 && INFO/DP<500 && \
    (INFO/FS<60 || INFO/FS=".") && INFO/MQ>40' input.vcf -o filtered.vcf

# Genotype-level filters
bcftools filter -i 'FMT/DP>10 && FMT/GQ>20' input.vcf -o filtered.vcf

# Remove filtered sites
bcftools view -f PASS input.vcf -o passed.vcf

# Keep only biallelic SNPs
bcftools view -m2 -M2 -v snps input.vcf -o biallelic_snps.vcf

# Check for missing values
bcftools filter -e 'QUAL="."' input.vcf.gz  # Exclude missing QUAL
bcftools filter -i 'INFO/DP!="."' input.vcf.gz  # Include only if DP exists
```

## bcftools view Filtering

**Goal:** Select variants by type, region, or sample using bcftools view.

**Approach:** Use type flags (-v/-V), region flags (-r/-R), and sample flags (-s/-S) for structured subsetting.

### Filter by Variant Type

```bash
# SNPs only
bcftools view -v snps input.vcf.gz -o snps.vcf.gz

# Indels only
bcftools view -v indels input.vcf.gz -o indels.vcf.gz

# Exclude SNPs
bcftools view -V snps input.vcf.gz -o no_snps.vcf.gz
```

### Filter by Region

```bash
bcftools view -r chr1:1000000-2000000 input.vcf.gz -o region.vcf.gz

# Multiple regions
bcftools view -r chr1:1000-2000,chr2:3000-4000 input.vcf.gz
```

### Filter by Samples

```bash
# Include samples
bcftools view -s sample1,sample2 input.vcf.gz -o subset.vcf.gz

# Exclude samples
bcftools view -s ^sample3,sample4 input.vcf.gz -o subset.vcf.gz
```

## Depth Filtering

**Goal:** Remove variants at extreme depth values that suggest mapping artifacts or duplications.

**Approach:** Calculate depth percentiles, then filter to the middle 90% of the distribution.

```bash
# Calculate depth percentiles
bcftools query -f '%DP\n' input.vcf | \
    sort -n | \
    awk '{a[NR]=$1} END {print "5th:", a[int(NR*0.05)], "95th:", a[int(NR*0.95)]}'

# Filter to middle 90% of depth distribution
bcftools filter -i 'INFO/DP>10 && INFO/DP<200' input.vcf -o depth_filtered.vcf
```

## Allele Frequency Filters

**Goal:** Filter variants by minor allele frequency and allelic balance.

**Approach:** Apply thresholds on INFO/AF for population frequency and AD ratios for heterozygote balance.

```bash
# Minor allele frequency filter (population data)
bcftools filter -i 'INFO/AF>0.01 && INFO/AF<0.99' input.vcf -o maf_filtered.vcf

# Allele balance for heterozygotes
bcftools filter -i 'GT="het" -> (AD[1]/(AD[0]+AD[1]) > 0.2 && AD[1]/(AD[0]+AD[1]) < 0.8)' \
    input.vcf -o ab_filtered.vcf
```

## Region-Based Filtering

**Goal:** Include or exclude variants based on genomic region annotations and artifact-prone regions.

**Approach:** Use bcftools view with BED files to exclude known problematic regions and restrict to targets. Always benchmark filtering stratified by genomic context.

Key BED resources for region filtering:
- **ENCODE exclusion list** (formerly blacklist): regions with anomalous signal across cell types (centromeres, telomeres, satellite repeats). Available at github.com/Boyle-Lab/Blacklist.
- **GIAB difficult regions**: stratified BED files for low-complexity, segmental duplications, tandem repeats, and other challenging contexts. Essential for honest benchmarking.
- **LCR (low-complexity region) BED files**: homopolymers, dinucleotide repeats, and other simple sequences where indel calling is unreliable. Heng Li's LCR-hs38 is widely used.

```bash
# Exclude ENCODE exclusion list regions
bcftools view -T ^ENCFF356LFX.bed input.vcf -o filtered.vcf

# Exclude low-complexity regions
bcftools view -T ^LCR-hs38.bed.gz input.vcf -o lcr_filtered.vcf

# Keep only exonic regions
bcftools view -R exons.bed input.vcf -o exonic.vcf

# Stratified evaluation against GIAB truth set
bcftools isec -p benchmark_dir filtered.vcf truth.vcf.gz
```

## Sample-Level Filtering

**Goal:** Remove low-quality samples and sites with excessive missing genotypes.

**Approach:** Compute per-sample missingness with bcftools stats, exclude failing samples, and filter sites by F_MISSING threshold.

```bash
# Missing genotype rate per sample
bcftools stats -s - input.vcf | grep ^PSC | cut -f3,14

# Filter samples with >10% missing
bcftools view -S good_samples.txt input.vcf -o sample_filtered.vcf

# Filter sites with >5% missing genotypes
bcftools filter -i 'F_MISSING<0.05' input.vcf -o site_filtered.vcf
```

## Variant Type-Specific Filters

**Goal:** Apply different quality thresholds to SNPs and indels.

**Approach:** Separate by type, apply type-appropriate thresholds (e.g., FS<60 for SNPs, FS<200 for indels), then merge back.

```bash
# Separate SNPs and indels for different filters
bcftools view -v snps input.vcf -o snps.vcf
bcftools view -v indels input.vcf -o indels.vcf

# Apply different thresholds
bcftools filter -i 'QUAL>30 && INFO/FS<60' snps.vcf -o snps_filtered.vcf
bcftools filter -i 'QUAL>30 && INFO/FS<200' indels.vcf -o indels_filtered.vcf

# Merge back
bcftools concat snps_filtered.vcf indels_filtered.vcf | bcftools sort -o merged.vcf
```

## Multi-Step Pipeline Filtering

**Goal:** Chain multiple filter criteria in a single pipeline with optional soft filter labels.

**Approach:** Pipe through successive bcftools filter commands, using -s for named soft filters and -f PASS for final hard extraction.

```bash
bcftools filter -e 'QUAL<30' input.vcf.gz | \
    bcftools filter -e 'INFO/DP<10' | \
    bcftools view -v snps -Oz -o filtered_snps.vcf.gz

# Named soft filters
bcftools filter -s 'LowQual' -e 'QUAL<30' input.vcf.gz | \
    bcftools filter -s 'LowDepth' -e 'INFO/DP<10' -Oz -o marked.vcf.gz

# Later, remove all filtered variants
bcftools view -f PASS marked.vcf.gz -Oz -o pass_only.vcf.gz
```

## Somatic Variant Filters

**Goal:** Filter somatic variants from tumor-normal calling with caller-specific criteria.

**Approach:** Use GATK FilterMutectCalls with contamination/segmentation tables, then apply additional bcftools thresholds on TLOD and VAF.

```bash
# Mutect2 specific
gatk FilterMutectCalls \
    -R reference.fa \
    -V mutect2_raw.vcf \
    --contamination-table contamination.table \
    --tumor-segmentation segments.table \
    -O mutect2_filtered.vcf

# Additional somatic filters
bcftools filter -i 'INFO/TLOD>6.3 && FMT/AF[0]>0.05 && FMT/DP[0]>20' \
    mutect2_filtered.vcf -o somatic_final.vcf
```

## Python Filtering (cyvcf2)

**Goal:** Filter variants programmatically in Python for custom multi-metric criteria.

**Approach:** Iterate with cyvcf2, evaluate QUAL/INFO/FORMAT fields per variant, and write passing records with Writer.

```python
from cyvcf2 import VCF, Writer
import numpy as np

vcf = VCF('input.vcf.gz')
writer = Writer('filtered.vcf', vcf)

for variant in vcf:
    qual = variant.QUAL or 0
    dp = variant.INFO.get('DP', 0)
    fs = variant.INFO.get('FS', 0)
    mq = variant.INFO.get('MQ', 0)

    if qual >= 30 and dp >= 10 and fs <= 60 and mq >= 40:
        writer.write_record(variant)

writer.close()
vcf.close()
```

### Filter by Genotype

```python
from cyvcf2 import VCF, Writer

vcf = VCF('input.vcf.gz')
writer = Writer('filtered.vcf', vcf)

for variant in vcf:
    # Keep if at least one sample is heterozygous
    # gt_types: 0=HOM_REF, 1=HET, 2=UNKNOWN, 3=HOM_ALT
    if 1 in variant.gt_types:
        writer.write_record(variant)

writer.close()
vcf.close()
```

### Filter by Sample Depth

```python
from cyvcf2 import VCF, Writer
import numpy as np

vcf = VCF('input.vcf.gz')
writer = Writer('filtered.vcf', vcf)

for variant in vcf:
    depths = variant.format('DP')
    if depths is not None and np.min(depths) >= 10:
        writer.write_record(variant)

writer.close()
vcf.close()
```

## Validate Filtering

**Goal:** Assess the impact of filtering on variant quality and retention using diagnostic metrics.

**Approach:** Compare before/after bcftools stats, check Ti/Tv and Het/Hom ratios against expected ranges, and verify known variant recovery.

Expected metric ranges for human samples:

| Metric | WGS | WES | Interpretation |
|--------|-----|-----|----------------|
| Ti/Tv | 2.0-2.1 | 2.8-3.3 | Below range suggests excess false positives; above range suggests over-filtering. WES is higher due to coding region enrichment (CpG transitions). |
| Het/Hom | 1.5-2.0 | 1.5-2.0 | Population-dependent. Outliers suggest contamination (high het) or inbreeding (low het). |
| Known variant % | >99% (dbSNP) | >99% | Low recovery of known variants indicates over-filtering. |

If Ti/Tv drops after filtering, the filters may be preferentially removing true transitions.

```bash
# Compare before/after stats
bcftools stats input.vcf > before_stats.txt
bcftools stats filtered.vcf > after_stats.txt

# Ti/Tv ratio
bcftools stats filtered.vcf | grep 'TSTV'

# Het/Hom ratio
bcftools stats -s - filtered.vcf | grep ^PSC | awk '{print $3, "het/hom:", $5/$6}'

# Count variants per filter label
bcftools query -f '%FILTER\n' filtered.vcf | sort | uniq -c

# Known variant recovery (requires dbSNP-annotated VCF)
bcftools view -i 'ID!="."' filtered.vcf | bcftools stats | grep 'number of SNPs'
```

## Common Filtering Pitfalls

- **Over-filtering rare variants**: Strict annotation thresholds disproportionately penalize rare variants, which have fewer supporting reads and therefore noisier annotations. Consider tiered filtering: relaxed thresholds for singletons, standard for common variants.
- **Applying SNP filters to indels**: SNPs and indels have fundamentally different annotation distributions (e.g., FS, SOR, ReadPosRankSum). Always separate by variant type before filtering.
- **Ignoring multiallelic sites**: Standard VQSR evaluates all alleles at a site together. For large cohorts, allele-specific VQSR (-AS flag) evaluates each allele independently, preventing a single poor allele from dragging down a good one.
- **Not verifying filter effects**: Always check Ti/Tv, Het/Hom, and known variant recovery rate before and after filtering. A filter that improves one metric while degrading another may be miscalibrated.
- **Filtering on depth alone**: DP filtering removes sites in collapsed segmental duplications that may harbor real variants. Combine DP with other annotations (MQ, MQRankSum) to distinguish true high-depth from mapping artifacts.
- **Not examining annotation distributions**: Plot histograms of each annotation (QD, FS, MQ) stratified by known true positives (e.g., HapMap sites) vs likely false positives before choosing thresholds. GATK defaults are population-level recommendations; specific datasets may warrant adjustment.

## Quick Reference

| Task | Command |
|------|---------|
| Quality filter | `bcftools filter -e 'QUAL<30' in.vcf.gz` |
| Depth filter | `bcftools filter -e 'INFO/DP<10' in.vcf.gz` |
| SNPs only | `bcftools view -v snps in.vcf.gz` |
| Indels only | `bcftools view -v indels in.vcf.gz` |
| PASS only | `bcftools view -f PASS in.vcf.gz` |
| Soft filter | `bcftools filter -s 'LowQ' -e 'QUAL<30' in.vcf.gz` |
| Region | `bcftools view -r chr1:1-1000 in.vcf.gz` |
| Biallelic | `bcftools view -m2 -M2 in.vcf.gz` |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `no such INFO tag` | Tag not in VCF | Check header with `bcftools view -h` |
| `syntax error` | Invalid expression | Check operator syntax (`\|\|` not `or`) |
| `empty output` | Filter too strict | Relax thresholds |

## Related Skills

- variant-calling/variant-calling - Variant calling with bcftools to generate VCF files
- variant-calling/gatk-variant-calling - GATK HaplotypeCaller and VQSR pipelines
- variant-calling/deepvariant - Deep learning variant calling as an alternative to GATK
- variant-calling/variant-annotation - Functional annotation after filtering
- variant-calling/vcf-statistics - Evaluate filter effects with concordance and summary metrics
- variant-calling/vcf-basics - VCF format fundamentals and field interpretation
- variant-calling/variant-normalization - Left-align and decompose before filtering for consistent comparisons
