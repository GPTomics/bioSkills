# Medaka Polishing and Variant Calling - Usage Guide

## Overview

Medaka uses neural networks to polish consensus sequences and call variants from Oxford Nanopore data. Models are trained on specific basecaller versions for optimal accuracy.

## When to Use This Skill

- You have a draft assembly from ONT data that needs polishing
- You want to call SNPs/indels from Nanopore reads
- You need to improve consensus accuracy for a genomic region

## Installation

```bash
# Conda (recommended)
conda create -n medaka -c conda-forge -c bioconda medaka
conda activate medaka

# Verify
medaka --version
```

## Model Selection

Critical: Use the model matching your basecaller.

```bash
# List available models
medaka tools list_models

# Common choices:
# - r1041_e82_400bps_sup_v5.0.0 (R10.4.1, SUP, latest)
# - r941_min_sup_g507 (R9.4.1, SUP)
```

## Basic Workflows

### Assembly Polishing

```bash
medaka_consensus \
    -i reads.fastq.gz \
    -d draft_assembly.fa \
    -o medaka_output \
    -t 4 \
    -m r1041_e82_400bps_sup_v5.0.0
```

### Variant Calling (Haploid)

```bash
# medaka v2.0+ uses medaka_variant for haploid samples
medaka_variant \
    -i reads.fastq.gz \
    -r reference.fa \
    -o output_dir \
    -m r1041_e82_400bps_sup_v5.0.0

# For diploid samples, use Clair3 instead (medaka diploid is deprecated)
```

## Common Issues

### Wrong Model

Symptoms: Poor polishing, unusual errors
Solution: Check basecaller version and select matching model

### Out of Memory

- Reduce batch size (-b 50)
- Process by region (--region chr1)
- Use fewer threads

### Slow Performance

- Use GPU if available
- Pre-align reads with minimap2
- Increase batch size if GPU memory allows

## Expected Accuracy

| Basecaller | Before Polish | After Polish |
|------------|---------------|--------------|
| SUP R10.4.1 | Q20 (~99%) | Q40+ (>99.99%) |
| HAC R10.4.1 | Q15-18 | Q30+ |

## Resources

- [Medaka GitHub](https://github.com/nanoporetech/medaka)
- [ONT Community](https://community.nanoporetech.com/)
