# Somatic Variant Pipeline Usage Guide

## Overview

Call somatic mutations from tumor-normal paired samples.

## Prerequisites

```bash
conda install -c bioconda gatk4 strelka bcftools
pip install pysam
```

## Pipeline Overview

```
Tumor BAM + Normal BAM
    ├── Preprocessing (MarkDuplicates, BQSR)
    ├── Variant Calling (Mutect2 or Strelka2)
    ├── Filtering (contamination, orientation bias)
    ├── Annotation (VEP, Funcotator)
    └── Final VCF
```

## Quick Start: Mutect2

```bash
# Call variants
gatk Mutect2 \
    -R reference.fa \
    -I tumor.bam \
    -I normal.bam \
    -normal normal_sample_name \
    --germline-resource gnomad.vcf.gz \
    --panel-of-normals pon.vcf.gz \
    -O somatic.vcf.gz

# Filter
gatk FilterMutectCalls \
    -R reference.fa \
    -V somatic.vcf.gz \
    -O filtered.vcf.gz
```

## Quick Start: Strelka2

```bash
# Configure
configureStrelkaSomaticWorkflow.py \
    --normalBam normal.bam \
    --tumorBam tumor.bam \
    --referenceFasta reference.fa \
    --runDir strelka_run

# Run
strelka_run/runWorkflow.py -m local -j 16
```

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

## Tips

- Always use matched normal when available
- Use gnomAD as germline resource
- Panel of normals reduces false positives
- Check tumor purity affects sensitivity

## Example Prompts

> "Call somatic mutations from my tumor-normal BAM pair using Mutect2"

> "Run the complete somatic variant pipeline with contamination estimation"

> "Create a panel of normals for somatic variant calling"

> "Annotate my somatic variants with VEP or Funcotator"

## See Also

- [GATK Mutect2 tutorial](https://gatk.broadinstitute.org/hc/en-us/articles/360035531132)
- [Strelka2 documentation](https://github.com/Illumina/strelka)
