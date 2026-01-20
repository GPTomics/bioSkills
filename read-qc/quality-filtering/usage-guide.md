# Quality Filtering - Usage Guide

## Overview
Quality filtering removes low-quality bases and reads that could cause downstream issues. Sliding window trimming preserves high-quality data while removing degraded 3' ends common in Illumina sequencing.

## Prerequisites
```bash
# Trimmomatic
conda install -c bioconda trimmomatic

# fastp
conda install -c bioconda fastp

# Cutadapt
pip install cutadapt
```

## Quick Start
Tell your AI agent what you want to do:
- "Filter low-quality reads from my FASTQ files"
- "Apply sliding window quality trimming to my data"
- "Remove poly-G artifacts from my NovaSeq data"

## Example Prompts

### Standard Filtering
> "Apply sliding window quality trimming with Q20 threshold to my paired-end reads"

> "Filter reads to keep only those with average quality above Q25 for variant calling"

### Platform-Specific
> "Remove poly-G tails from my NovaSeq data using fastp"

> "My NextSeq data has quality issues, apply appropriate filtering"

### Parameter Tuning
> "Too many reads are being discarded, relax the quality filtering parameters"

> "I need high-accuracy reads for variant calling, use stringent quality filtering"

## What the Agent Will Do
1. Assess current quality distribution from FastQC reports
2. Select appropriate tool and parameters based on application
3. Apply sliding window or global quality filtering
4. Set minimum read length threshold
5. Generate filtered output with quality metrics

## Choosing Parameters

### Quality Threshold by Application

| Application | Threshold | Notes |
|-------------|-----------|-------|
| General | Q20 | 1% error rate |
| Variant calling | Q25-30 | Need high accuracy |
| Assembly | Q15-20 | Quantity matters |
| RNA-seq | Q20 | Standard |

### Minimum Length by Read Type

| Read Type | Min Length | Notes |
|-----------|------------|-------|
| Short reads (150bp) | 36-50 | ~25-35% of original |
| Short reads (100bp) | 30-36 | |
| Long inserts | 50+ | May need longer |

## Tips
- Use `--cut_right` in fastp for sliding window behavior similar to Trimmomatic
- NovaSeq/NextSeq require poly-G trimming due to two-color chemistry
- Start with standard parameters and adjust based on read loss
- Higher thresholds for variant calling, lower for assembly where quantity matters
- Always check FastQC before and after filtering to verify improvement

## Resources
- [Trimmomatic Manual](http://www.usadellab.org/cms/?page=trimmomatic)
- [fastp Documentation](https://github.com/OpenGene/fastp)
- [Quality Score Tutorial](https://www.drive5.com/usearch/manual/quality_score.html)
