# Variant Filtering Best Practices - Usage Guide

## Overview

Variant filtering decides which errors a callset keeps, not whether it has them. Site-level filtering (VQSR, VETS, NVScoreVariants, or hard thresholds) asks whether a variant SITE is real; genotype-level filtering (GQ/DP/allele balance) asks whether an individual sample's genotype at a passing site is trustworthy. The two are orthogonal and both are needed, site first. VQSR and hard filters fail in opposite regimes -- VQSR collapses on small/exome cohorts, hard filters discard real variants at scale -- so the right method depends on cohort size, sequencing platform, and organism. None of the common mistakes throw an error; the VCF stays structurally valid while the numbers are silently wrong.

## Prerequisites

```bash
# bcftools
conda install -c bioconda bcftools

# GATK
conda install -c bioconda gatk4

# cyvcf2 (Python filtering)
pip install cyvcf2
```

## Quick Start

Tell your AI agent what you want to do:
- "Apply GATK hard filters to my SNP and indel calls separately"
- "Filter variants with bcftools for quality and depth"
- "Set up VQSR for my whole genome cohort"
- "Decide whether to use VQSR or hard filters for my exome"
- "Null out low-GQ genotypes before I compute cohort HWE"
- "Filter somatic variants from my tumor-normal calls"

## Example Prompts

### Choosing a Method
> "I have a single exome -- should I run VQSR or hard filters?"

> "My cohort has 200 jointly-genotyped WGS samples; set up allele-specific VQSR"

> "CNNScoreVariants is deprecated in my GATK build -- what should I use for one sample?"

### Hard Filtering
> "Apply GATK hard filters to my SNP calls with the recommended thresholds"

> "Apply standard hard filter thresholds for indels and explain why they differ from SNPs"

> "Reproduce the GATK hard filter in bcftools so hom-alt sites don't get dropped"

### VQSR
> "Run VQSR on my whole genome cohort and apply the 99.7 truth-sensitivity tranche"

> "Explain what VQSLOD and truth-sensitivity tranches mean in my recalibrated VCF"

### Genotype and Population Filtering
> "Set low-confidence genotypes to no-call before I compute missingness and HWE"

> "Filter my cohort VCF by MAF and per-sample missingness"

> "Remove variants in ENCODE exclusion-list and low-complexity regions"

### Somatic Filtering
> "Set up filtering for somatic variant calls from Mutect2"

### Understanding Metrics
> "Explain what QD, FS, SOR, and MQRankSum mean in my VCF"

> "My WES Ti/Tv is 2.1 after filtering -- is that a problem?"

## What the Agent Will Do

1. Assess dataset characteristics (cohort size, variant type, caller, organism, WGS vs WES)
2. Select a site-level method: hard filters, VQSR, allele-specific VQSR, VETS, or NVScoreVariants (or none, for DeepVariant/DRAGEN output)
3. Split SNPs and indels and apply type-appropriate thresholds, then merge
4. Guard RankSum annotations so hom-alt sites are not silently removed
5. Apply genotype-level GQ/DP/allele-balance no-calls before cohort QC
6. Exclude artifact-prone regions if BED files are provided
7. Validate with Ti/Tv, Het/Hom, and known-variant recovery

## Relocated Reference Recipes

The core decision content, VQSR/hard-filter commands, and the hom-alt trap live in SKILL.md. The subsetting and cohort-QC recipes below are the fuller reference set.

### Subset by Type, Region, or Sample

```bash
bcftools view -v snps input.vcf.gz -o snps.vcf.gz        # SNPs only
bcftools view -v indels input.vcf.gz -o indels.vcf.gz    # indels only
bcftools view -V snps input.vcf.gz -o no_snps.vcf.gz     # exclude SNPs
bcftools view -m2 -M2 -v snps input.vcf.gz -o biallelic_snps.vcf.gz  # biallelic SNPs

bcftools view -r chr1:1000000-2000000 input.vcf.gz -o region.vcf.gz  # region
bcftools view -s sample1,sample2 input.vcf.gz -o subset.vcf.gz       # keep samples
bcftools view -s ^sample3 input.vcf.gz -o subset.vcf.gz              # drop samples
```

### Depth-Distribution Filtering

```bash
# Inspect the depth distribution before choosing cutoffs
bcftools query -f '%DP\n' input.vcf | sort -n | \
    awk '{a[NR]=$1} END {print "5th:", a[int(NR*0.05)], "95th:", a[int(NR*0.95)]}'

# Filter to the middle of the depth distribution (extreme depth => collapsed repeats/CNV)
bcftools filter -i 'INFO/DP>10 && INFO/DP<200' input.vcf -o depth_filtered.vcf
```

### Allele-Frequency and Allele-Balance Filters

```bash
# Minor allele frequency (population data)
bcftools filter -i 'INFO/AF>0.01 && INFO/AF<0.99' input.vcf -o maf_filtered.vcf

# Allele balance at hets: a true het is ~0.5 alt fraction; far from 0.5 is suspect
bcftools filter -i 'GT="het" && (AD[1]/(AD[0]+AD[1]) > 0.2 && AD[1]/(AD[0]+AD[1]) < 0.8)' \
    input.vcf -o ab_filtered.vcf
```

### Per-Sample and Per-Site Missingness

```bash
# Per-sample missingness (drop bad samples, then recompute per-variant call rate)
bcftools stats -s - input.vcf | grep ^PSC | cut -f3,14
bcftools view -S good_samples.txt input.vcf -o sample_filtered.vcf

# Drop sites with >5% missing genotypes
bcftools filter -i 'F_MISSING<0.05' input.vcf -o site_filtered.vcf
```

Variant-level and sample-level missingness are coupled: drop the worst samples, recompute per-variant call rate, drop the worst variants, iterate -- do not apply both thresholds in one pass. Remember `./.` (no-call) is not `0/0` (hom-ref); only no-calls count as missing.

### Multi-Step Soft-Filter Pipeline

```bash
# Label failures at each stage, then extract PASS at the end
bcftools filter -s 'LowQual' -e 'QUAL<30' input.vcf.gz | \
    bcftools filter -s 'LowDepth' -e 'INFO/DP<10' -Oz -o marked.vcf.gz
bcftools view -f PASS marked.vcf.gz -Oz -o pass_only.vcf.gz
```

### Region Exclusion BED Resources

```bash
bcftools view -T ^ENCFF356LFX.bed input.vcf -o filtered.vcf   # ENCODE exclusion list
bcftools view -T ^LCR-hs38.bed.gz input.vcf -o lcr_filtered.vcf  # low-complexity regions
bcftools isec -p benchmark_dir filtered.vcf.gz truth.vcf.gz    # stratified GIAB benchmark
```

## Tips

- Always separate SNPs and indels before filtering -- their annotation distributions differ.
- Guard every RankSum term with `|| INFO/X = "."` in bcftools, or hom-alt sites silently vanish (RankSum is undefined at hom-alt sites).
- Do not run VQSR on a single exome or panel (too few variants; a single deep WGS is fine); the Gaussian mixture is non-identifiable there. Use hard filters, VETS, or NVScoreVariants.
- Never use DP as a VQSR annotation for exomes -- capture depth tracks bait design, not truth.
- Do not re-apply GATK hard filters on DeepVariant/DRAGEN output; filter on the caller's own calibrated fields.
- Ti/Tv ~2.0-2.1 (WGS) or ~3.0-3.3 (WES); below range signals excess false positives.
- Apply genotype-level no-calls before computing cohort missingness/HWE, or garbage genotypes drive spurious deviation.
- Do not over-filter rare variants: a real singleton looks statistically like an artifact, so aggressive site filters preferentially delete true rare alleles. For rare-variant work apply tiered filtering (looser thresholds with an orthogonal check such as manual IGV review or replication) rather than one blanket cutoff.
- cyvcf2 genotype codes for `variant.gt_types`: 0=HOM_REF, 1=HET, 2=UNKNOWN (missing), 3=HOM_ALT.

## Filter Method Selection

| Dataset | Recommended approach |
|---------|---------------------|
| Single exome / gene panel | Hard filters, VETS, or NVScoreVariants (not VQSR) |
| Single deep WGS, or ~30+ jointly-genotyped exomes (human) | VQSR |
| Biobank-scale cohort | Allele-specific VQSR (`-AS`) or VETS |
| Non-model organism (no truth sets) | Hard filters only |
| DeepVariant / DRAGEN output | Filter on the caller's own fields |
| Somatic / tumor-normal | FilterMutectCalls + caller-specific thresholds |

## Related Skills

- variant-calling/gatk-variant-calling - HaplotypeCaller and joint genotyping upstream of VQSR
- variant-calling/deepvariant - Deep-learning caller whose output needs no separate site filter
- variant-calling/variant-normalization - Left-align and decompose before filtering
- variant-calling/vcf-statistics - Ti/Tv, het/hom, contamination, and relatedness QC
- variant-calling/vcf-basics - VCF field interpretation and PL/GQ/QUAL relationships
