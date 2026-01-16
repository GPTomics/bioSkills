# Quality Filtering - Usage Guide

## Overview

Quality filtering removes low-quality bases and reads that could cause downstream issues. Sliding window trimming preserves high-quality data while removing degraded 3' ends common in Illumina sequencing.

## When to Use This Skill

- FastQC shows quality drop at read ends
- Preparing reads for variant calling (requires high accuracy)
- Removing low-quality reads before assembly
- Processing NovaSeq/NextSeq data with poly-G artifacts

## Installation

```bash
# Trimmomatic
conda install -c bioconda trimmomatic

# fastp
conda install -c bioconda fastp

# Cutadapt
pip install cutadapt
```

## Common Workflows

### Standard Quality Trimming (Trimmomatic)

```bash
trimmomatic PE -threads 4 \
    R1.fq.gz R2.fq.gz \
    R1_clean.fq.gz R1_unpaired.fq.gz \
    R2_clean.fq.gz R2_unpaired.fq.gz \
    SLIDINGWINDOW:4:20 MINLEN:36
```

### Modern Workflow (fastp)

```bash
fastp -i R1.fq.gz -I R2.fq.gz \
      -o R1_clean.fq.gz -O R2_clean.fq.gz \
      --cut_right -q 20 -l 36
```

## Choosing Parameters

### Window Size

- **4 bp**: Standard, balances sensitivity and specificity
- **5-10 bp**: More conservative, smoother signal

### Quality Threshold

| Application | Threshold | Notes |
|-------------|-----------|-------|
| General | Q20 | 1% error rate |
| Variant calling | Q25-30 | Need high accuracy |
| Assembly | Q15-20 | Quantity matters |
| RNA-seq | Q20 | Standard |

### Minimum Length

| Read Type | Min Length | Notes |
|-----------|------------|-------|
| Short reads (150bp) | 36-50 | ~25-35% of original |
| Short reads (100bp) | 30-36 | |
| Long inserts | 50+ | May need longer |

## Troubleshooting

### Too Many Reads Lost

Relax parameters:
```bash
# Lower quality threshold
SLIDINGWINDOW:4:15

# Shorter minimum length
MINLEN:25
```

### Quality Still Poor After Filtering

Increase stringency:
```bash
# Higher threshold
SLIDINGWINDOW:4:25

# Trim both ends
LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20
```

### Poly-G Artifacts (NovaSeq)

Use fastp with poly-G trimming:
```bash
fastp -i in.fq -o out.fq --trim_poly_g
```

## Resources

- [Trimmomatic Manual](http://www.usadellab.org/cms/?page=trimmomatic)
- [fastp Documentation](https://github.com/OpenGene/fastp)
- [Quality Score Tutorial](https://www.drive5.com/usearch/manual/quality_score.html)
