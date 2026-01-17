# GATK Variant Calling Usage Guide

This guide covers germline variant calling using GATK HaplotypeCaller.

## When to Use GATK

- Production-quality variant calls
- Cohort analysis with joint genotyping
- When VQSR is possible (many variants)
- Following GATK Best Practices

## Workflow Options

### Single Sample
1. Mark duplicates (Picard/samtools)
2. BQSR (optional but recommended)
3. HaplotypeCaller
4. Hard filtering

### Cohort (Recommended)
1. Mark duplicates per sample
2. BQSR per sample
3. HaplotypeCaller with -ERC GVCF per sample
4. GenomicsDBImport or CombineGVCFs
5. GenotypeGVCFs
6. VQSR or hard filtering

## GVCF vs VCF

- **VCF**: Final variant calls, reference sites not included
- **GVCF**: Includes reference blocks for joint genotyping
- Always use GVCF for cohorts to capture reference confidence

## VQSR vs Hard Filtering

### Use VQSR when:
- Whole genome sequencing
- Large cohort (many variants)
- Resource files available

### Use Hard Filtering when:
- Exome/panel sequencing
- Single sample or small cohort
- Non-model organism
- VQSR fails (not enough variants)

## Tool Requirements

```bash
# GATK 4.x
conda install -c bioconda gatk4

# Or download from Broad
# https://github.com/broadinstitute/gatk/releases
```

## Resource Files

Download from GATK Resource Bundle:
- gs://genomics-public-data/resources/broad/hg38/v0/

Key files:
- Homo_sapiens_assembly38.fasta
- Homo_sapiens_assembly38.dbsnp138.vcf
- hapmap_3.3.hg38.vcf.gz
- 1000G_omni2.5.hg38.vcf.gz
- Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
