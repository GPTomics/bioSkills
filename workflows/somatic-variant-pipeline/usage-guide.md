# Somatic Variant Pipeline Usage Guide

## Overview

Call somatic mutations from tumor-normal paired samples using Mutect2 or Strelka2, with contamination estimation, orientation bias filtering, and annotation.

## Prerequisites

```bash
conda install -c bioconda gatk4 strelka bcftools
pip install pysam
```

## Quick Start

Tell your AI agent what you want to do:
- "Call somatic mutations from my tumor-normal BAM pair"
- "Run Mutect2 with contamination estimation and orientation bias filtering"
- "Create a panel of normals from my normal samples"
- "Annotate somatic variants with VEP"

## Complete Mutect2 Workflow

### 1. Calculate Contamination

```bash
gatk GetPileupSummaries \
    -I tumor.bam \
    -V gnomad.vcf.gz \
    -L intervals.bed \
    -O tumor_pileups.table

gatk CalculateContamination \
    -I tumor_pileups.table \
    -O contamination.table
```

### 2. Learn Read Orientation Model

```bash
gatk LearnReadOrientationModel \
    -I f1r2.tar.gz \
    -O read-orientation-model.tar.gz
```

### 3. Filter with Contamination

```bash
gatk FilterMutectCalls \
    -R reference.fa \
    -V somatic.vcf.gz \
    --contamination-table contamination.table \
    --ob-priors read-orientation-model.tar.gz \
    -O filtered.vcf.gz
```

## Annotation

```bash
# VEP
vep -i filtered.vcf.gz -o annotated.vcf \
    --cache --assembly GRCh38 \
    --vcf --symbol --everything

# Funcotator
gatk Funcotator \
    -R reference.fa \
    -V filtered.vcf.gz \
    -O funcotated.vcf \
    --data-sources-path funcotator_dataSources \
    --output-file-format VCF
```

## Tumor-Only Mode

When normal sample unavailable:

```bash
gatk Mutect2 \
    -R reference.fa \
    -I tumor.bam \
    --germline-resource gnomad.vcf.gz \
    --panel-of-normals pon.vcf.gz \
    -O tumor_only.vcf.gz
```

## Panel of Normals

```bash
# Create PON from normal samples
gatk CreateSomaticPanelOfNormals \
    -V normal1.vcf.gz \
    -V normal2.vcf.gz \
    -V normal3.vcf.gz \
    -O pon.vcf.gz
```

## Quality Metrics

```bash
# Variant statistics
bcftools stats filtered.vcf.gz > stats.txt

# Count by type
bcftools view -f PASS filtered.vcf.gz | bcftools stats -
```

## Example Prompts

> "Call somatic mutations from my tumor-normal BAM pair using Mutect2"

> "Run the complete somatic variant pipeline with contamination estimation"

> "Create a panel of normals for somatic variant calling"

> "Annotate my somatic variants with VEP or Funcotator"

## What the Agent Will Do

1. Run Mutect2 or Strelka2 on tumor-normal paired BAMs
2. Estimate cross-sample contamination and learn orientation bias model
3. Apply FilterMutectCalls with contamination and orientation bias corrections
4. Extract PASS somatic variants
5. Annotate with VEP or Funcotator for functional interpretation

## Tips

- Always use matched normal when available -- tumor-only mode has higher false positive rate
- Use gnomAD as germline resource to filter common germline variants
- Panel of normals (40+ normals from same platform) reduces systematic artifacts
- Tumor purity affects sensitivity -- low-purity samples need higher sequencing depth
- Consensus calling with 2+ callers (Mutect2 + Strelka2) improves accuracy

## Related Skills

- variant-calling/gatk-variant-calling - Germline variant calling
- variant-calling/filtering-best-practices - Filtering strategies
- variant-calling/variant-annotation - VEP/SnpEff annotation
- variant-calling/structural-variant-calling - Somatic SV detection with Manta
- copy-number/cnvkit-analysis - Somatic CNV calling
