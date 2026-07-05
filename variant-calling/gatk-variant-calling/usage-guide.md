# GATK Variant Calling - Usage Guide

## Overview

Germline SNP and indel calling with GATK HaplotypeCaller, the reference implementation for production variant calling. HaplotypeCaller does not trust the input alignment: in any region showing signal it locally reassembles candidate haplotypes and realigns reads against them with a PairHMM, which is why it beats pileup callers (bcftools) on indels and complex loci. For cohorts it emits per-sample GVCFs (`-ERC GVCF`) that are jointly genotyped later, decoupling expensive per-sample discovery from cheap cohort genotyping.

## Prerequisites

```bash
# GATK 4.x
conda install -c bioconda gatk4

# Or download the release jar from Broad
# https://github.com/broadinstitute/gatk/releases
```

Inputs: a coordinate-sorted, indexed BAM with read groups and duplicates marked (Picard/GATK MarkDuplicates or samtools markdup); a reference FASTA with `.fai` and `.dict`; for BQSR/VQSR, the GATK resource-bundle known-sites and truth files.

## Quick Start

Tell your AI agent what you want to do:
- "Call germline variants from my aligned BAM with HaplotypeCaller"
- "Run HaplotypeCaller in GVCF mode so I can joint-genotype the cohort later"
- "Use DRAGEN-GATK mode for single-sample calling without BQSR"
- "Decide whether I still need BQSR on my NovaSeq data"
- "Call variants on chrX correctly for a male sample"
- "Call mitochondrial heteroplasmy from my BAM"
- "Check my sample for contamination before I trust the calls"

## Example Prompts

### Single Sample Calling
> "Call variants from sample.bam using GATK HaplotypeCaller"

> "Run HaplotypeCaller on my exome with the capture-regions interval list"

### GVCF for Joint Calling
> "Generate a GVCF from each BAM in my cohort so I can joint-genotype them"

> "Explain why I should call per-sample GVCFs instead of calling the whole cohort together"

### DRAGEN and BQSR decisions
> "Run HaplotypeCaller in DRAGEN mode and hard-filter on QUAL"

> "Is BQSR still worth running on my binned-quality NovaSeq data?"

### Edge cases
> "Set the right ploidy to call non-PAR chrX and chrY in a 46,XY sample"

> "Call mitochondrial variants with Mutect2 mitochondria mode and the shifted reference"

> "Estimate cross-sample contamination with VerifyBamID2 or CHARR before variant calling"

### Filtering
> "Apply VQSR to my whole-genome cohort calls" (see filtering-best-practices)

> "Hard-filter my small exome callset with GATK-recommended thresholds"

## What the Agent Will Do

1. Verify the input BAM is sorted, indexed, has read groups, and has duplicates marked; verify the reference has `.fai`/`.dict`.
2. Decide the mode from context: standard vs `--dragen-mode`, direct VCF vs `-ERC GVCF` for cohorts, and the correct `--sample-ploidy` for non-diploid/sex-chromosome/mitochondrial cases.
3. Optionally run BQSR (BaseRecalibrator + ApplyBQSR) for standard mode -- skipped in DRAGEN mode and often low-impact on binned-quality data.
4. Run HaplotypeCaller (or Mutect2 for mitochondria/somatic), scattering by interval for speed.
5. Genotype (GenotypeGVCFs) for single samples; hand cohorts to the joint-calling workflow.
6. Choose a filtering method (VQSR/VETS, AS_VQSR, hard filters, or DRAGEN QUAL) and apply it, then compute Ti/Tv and callset statistics.

## Relocated recipes

The SKILL.md keeps the decision content lean; the fuller command recipes it points to are collected here.

### Hard filtering (single sample / small cohort / non-model organism)

SNPs and indels are filtered separately because their annotation distributions differ. These are GATK-recommended starting points, not universal truth -- plot annotation histograms for known-true/false sites and adjust.

```bash
gatk SelectVariants -R reference.fa -V cohort.vcf.gz --select-type-to-include SNP -O snps.vcf.gz
gatk SelectVariants -R reference.fa -V cohort.vcf.gz --select-type-to-include INDEL -O indels.vcf.gz

gatk VariantFiltration -R reference.fa -V snps.vcf.gz -O snps.filtered.vcf.gz \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "FS > 60.0" --filter-name "FS60" \
    --filter-expression "MQ < 40.0" --filter-name "MQ40" \
    --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    --filter-expression "SOR > 3.0" --filter-name "SOR3"

gatk VariantFiltration -R reference.fa -V indels.vcf.gz -O indels.filtered.vcf.gz \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "FS > 200.0" --filter-name "FS200" \
    --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    --filter-expression "SOR > 10.0" --filter-name "SOR10"

gatk MergeVcfs -I snps.filtered.vcf.gz -I indels.filtered.vcf.gz -O cohort.filtered.vcf.gz
```

Threshold meanings: QD (variant quality per unit depth), FS/SOR (strand bias; indels tolerate more, hence higher thresholds), MQ (mapping ambiguity in paralogous/segdup regions), MQRankSum (alt reads map worse than ref), ReadPosRankSum (variant clusters at read ends). Full threshold rationale lives in filtering-best-practices.

### VQSR (large human WGS cohort)

VQSR fits a Gaussian mixture model over site annotations using truth/training resources; run SNP and INDEL modes separately. It needs many variants overlapping the resources to converge -- use hard filtering otherwise.

```bash
gatk VariantRecalibrator -R reference.fa -V cohort.vcf.gz \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap.vcf.gz \
    --resource:omni,known=false,training=true,truth=false,prior=12.0 omni.vcf.gz \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp.vcf.gz \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode SNP -O snp.recal --tranches-file snp.tranches

gatk ApplyVQSR -R reference.fa -V cohort.vcf.gz -O cohort.snp_recal.vcf.gz \
    --recal-file snp.recal --tranches-file snp.tranches \
    --truth-sensitivity-filter-level 99.5 -mode SNP
```

Add `-AS` (allele-specific VQSR) for large cohorts with many multiallelic sites; requires AS_ annotations (`-G AS_StandardAnnotation`). Indel mode uses the Mills gold-standard resource and commonly `--max-gaussians 4`. GATK is deprecating VQSR toward VETS (ExtractVariantAnnotations -> TrainVariantAnnotationsModel -> ScoreVariantAnnotations) -- see filtering-best-practices.

## Tips

- Always use `-ERC GVCF` when joint genotyping is possible later; the `<NON_REF>` allele is what lets a variant found in another sample be evaluated in this one.
- Mark duplicates before calling in every mode. In DRAGEN mode skip BQSR entirely.
- For WGS, expect a Ti/Tv ratio around 2.0-2.1 for high-quality SNP calls (lower suggests false positives).
- Set `--sample-ploidy` explicitly for pooled, polyploid, and non-PAR sex-chromosome calling; the diploid default emits impossible genotypes otherwise.
- Use Mutect2 `--mitochondria-mode` (not HaplotypeCaller) for mtDNA; heteroplasmy is continuous-VAF like subclonal somatic variation.
- Gate the callset on a contamination estimate (VerifyBamID2 or CHARR); even 1-3% contamination fabricates false hets.
- The PairHMM dominates runtime; scatter by interval and use `--native-pair-hmm-threads`.

## Resource files

Download from the GATK Resource Bundle (`gs://gcp-public-data--broad-references/hg38/v0/`):
- Homo_sapiens_assembly38.fasta (+ .fai, .dict)
- Homo_sapiens_assembly38.dbsnp138.vcf (known sites; prior 2.0)
- hapmap_3.3.hg38.vcf.gz (training/truth SNPs; prior 15.0)
- 1000G_omni2.5.hg38.vcf.gz (training SNPs; prior 12.0)
- 1000G_phase1.snps.high_confidence.hg38.vcf.gz (training SNPs; prior 10.0)
- Mills_and_1000G_gold_standard.indels.hg38.vcf.gz (training/truth indels; prior 12.0)

## Related Skills

- variant-calling/joint-calling - Cohort consolidation and joint genotyping at scale
- variant-calling/filtering-best-practices - Full VQSR/VETS and hard-filter recipes with rationale
- variant-calling/deepvariant - CNN-based alternative caller
- variant-calling/variant-calling - Lightweight bcftools mpileup/call
- variant-calling/variant-normalization - Normalize indels after calling
- variant-calling/vcf-basics - View and query resulting VCF files
- read-alignment/bwa-alignment - Align reads before variant calling
