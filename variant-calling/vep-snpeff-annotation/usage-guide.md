# VEP/SnpEff Annotation - Usage Guide

## Overview

VEP (Ensembl Variant Effect Predictor), SnpEff, and ANNOVAR provide comprehensive variant annotation beyond bcftools csq. They predict functional consequences, add pathogenicity scores, population frequencies, and clinical significance.

## When to Use This Skill

- Clinical variant interpretation
- Filtering to pathogenic/likely pathogenic variants
- Adding SIFT/PolyPhen/CADD predictions
- Population frequency filtering (gnomAD)
- HGVS nomenclature for reporting
- Research requiring detailed consequence annotations

## Tool Selection

| Use Case | Recommended Tool |
|----------|------------------|
| Clinical/research | VEP (Ensembl integration, plugins) |
| Quick annotation | SnpEff (fastest) |
| Custom databases | ANNOVAR (flexible) |
| Production pipeline | SnpEff (speed) or VEP (depth) |

## Installation

```bash
# VEP
conda install -c bioconda ensembl-vep
vep_install -a cf -s homo_sapiens -y GRCh38

# SnpEff
conda install -c bioconda snpeff
snpEff download GRCh38.105

# ANNOVAR (requires registration)
# Download from https://annovar.openbioinformatics.org/
```

## Quick Start

### VEP

```bash
vep -i input.vcf -o output.vcf --vcf --cache --offline --everything --pick
```

### SnpEff

```bash
snpEff ann GRCh38.105 input.vcf > output.vcf
```

### ANNOVAR

```bash
table_annovar.pl input.vcf humandb/ -buildver hg38 \
    -protocol refGene,gnomad30_genome,clinvar -operation g,f,f -vcfinput
```

## Key Concepts

### Impact Levels

| Level | Examples | Clinical Relevance |
|-------|----------|-------------------|
| HIGH | Frameshift, stop gained, splice | Likely functional |
| MODERATE | Missense, inframe indel | Possible functional |
| LOW | Synonymous, splice region | Unlikely functional |
| MODIFIER | Intron, UTR, intergenic | Usually benign |

### Pathogenicity Thresholds

| Score | Deleterious | Benign |
|-------|-------------|--------|
| SIFT | < 0.05 | >= 0.05 |
| PolyPhen | > 0.85 | < 0.15 |
| CADD | > 20 | < 10 |

## Common Workflows

### Clinical Filtering

```bash
# VEP with clinical databases
vep -i input.vcf -o output.vcf --vcf --cache --offline \
    --everything --pick \
    --plugin CADD,cadd.tsv.gz

# Filter to variants of interest
bcftools view -i 'INFO/CSQ~"HIGH" || INFO/CSQ~"MODERATE"' output.vcf
```

### Research Annotation

```bash
# Full SnpEff + SnpSift pipeline
snpEff ann GRCh38.105 input.vcf | \
    SnpSift annotate dbsnp.vcf.gz | \
    SnpSift annotate gnomad.vcf.gz | \
    SnpSift annotate clinvar.vcf.gz > annotated.vcf
```

## Resources

- [VEP Documentation](https://www.ensembl.org/info/docs/tools/vep/index.html)
- [SnpEff Documentation](https://pcingola.github.io/SnpEff/)
- [ANNOVAR Documentation](https://annovar.openbioinformatics.org/)
