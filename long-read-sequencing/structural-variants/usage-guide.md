# Structural Variant Detection - Usage Guide

## Overview

Long reads excel at detecting structural variants (deletions, insertions, inversions, duplications) that are difficult or impossible to detect with short reads. Multiple callers exist with different strengths.

## Tool Comparison

| Tool | Best For | Notes |
|------|----------|-------|
| Sniffles2 | General purpose | Most popular, population calling |
| cuteSV | High sensitivity | Platform-specific settings |
| SVIM | ONT data | Good insertion detection |
| pbsv | PacBio HiFi | PacBio-optimized |

## When to Use This Skill

- You have long-read alignments and want to find SVs
- You're looking for insertions/deletions >50bp
- You need to detect inversions or complex rearrangements
- You're doing disease research involving SVs

## Installation

```bash
conda install -c bioconda sniffles cutesv svim
```

## Basic Workflow

```bash
# 1. Align reads
minimap2 -ax map-ont ref.fa reads.fq.gz | samtools sort -o aligned.bam
samtools index aligned.bam

# 2. Call SVs
sniffles --input aligned.bam --vcf svs.vcf --reference ref.fa

# 3. Filter
bcftools filter -i 'QUAL>=20' svs.vcf > svs.filtered.vcf
```

## Best Practices

### Use Reference for Insertions

Sniffles2 needs reference to output insertion sequences:
```bash
sniffles --input aligned.bam --vcf svs.vcf --reference ref.fa
```

### Platform-Specific Settings

- ONT: Higher error rate, adjust parameters
- HiFi: Lower error, can use stricter thresholds

### Tandem Repeat Annotation

Improves accuracy in repetitive regions:
```bash
sniffles --input aligned.bam --vcf svs.vcf --tandem-repeats tandem_repeats.bed
```

## Common Issues

### Too Few SVs

- Check coverage (need >10x)
- Lower min support threshold
- Verify alignment quality

### Too Many False Positives

- Increase quality threshold
- Increase min support
- Filter repetitive regions

## Resources

- [Sniffles2 GitHub](https://github.com/fritzsedlazeck/Sniffles)
- [cuteSV GitHub](https://github.com/tjiangHIT/cuteSV)
