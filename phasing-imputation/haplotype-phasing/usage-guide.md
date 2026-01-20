# Haplotype Phasing - Usage Guide

## Overview
Haplotype phasing determines which alleles are inherited together on each chromosome, essential for genotype imputation, haplotype-based association tests, and compound heterozygosity detection.

## Prerequisites
```bash
# Beagle (download jar)
wget https://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar

# SHAPEIT5 (for large datasets)
conda install -c bioconda shapeit5

# bcftools for data preparation
conda install -c bioconda bcftools

# Genetic maps (required for phasing)
wget https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip
unzip plink.GRCh38.map.zip -d genetic_maps/
```

## Quick Start
Tell your AI agent what you want to do:
- "Phase my VCF file using Beagle"
- "Set up a phasing pipeline for all chromosomes"
- "Phase my large biobank data with SHAPEIT5"

## Example Prompts
### Basic Phasing
> "Phase my genotype data in study.vcf.gz using Beagle with GRCh38 genetic maps"

### Large Dataset Phasing
> "Use SHAPEIT5 two-stage approach to phase my 100,000 sample dataset"

### Reference-Assisted Phasing
> "Phase my study data using 1000 Genomes as a reference panel for better rare variant phasing"

### Quality Evaluation
> "Evaluate phasing quality using trio data to calculate switch error rates"

## What the Agent Will Do
1. Prepare input data by filtering to biallelic SNPs and checking for strand issues
2. Download required genetic maps for the appropriate genome build
3. Run phasing tool (Beagle or SHAPEIT5) per chromosome
4. Merge phased chromosomes into final output
5. Evaluate phasing quality if trio data is available

## Tips
- Always use genetic maps - they significantly improve phasing accuracy
- Process chromosomes in parallel since each is independent
- For datasets >50,000 samples, use SHAPEIT5's two-stage approach
- Filter to biallelic SNPs before phasing (Beagle handles these best)
- Using a reference panel improves accuracy for rare variants
- Memory errors can be resolved by increasing heap size or phasing in smaller chunks

## Choosing a Phasing Tool

| Tool | Best For | Memory | Speed |
|------|----------|--------|-------|
| Beagle 5.4 | General purpose, <50k samples | Medium | Fast |
| SHAPEIT5 | Large biobanks (>50k samples) | Low | Very fast |
| Eagle2 | Phasing with ref panel | Medium | Fast |

## Beagle Workflow

```bash
# Prepare input
bcftools view -m2 -M2 -v snps input.vcf.gz -Oz -o biallelic.vcf.gz
bcftools index biallelic.vcf.gz

# Phase per chromosome
for chr in {1..22}; do
    java -Xmx16g -jar beagle.jar \
        gt=biallelic.vcf.gz \
        chrom=chr${chr} \
        map=genetic_maps/plink.chr${chr}.GRCh38.map \
        out=phased_chr${chr} \
        nthreads=8
done

# Merge chromosomes
ls phased_chr*.vcf.gz > filelist.txt
bcftools concat -f filelist.txt -Oz -o phased.vcf.gz
bcftools index phased.vcf.gz
```

## SHAPEIT5 Workflow (Large Datasets)

```bash
# For datasets >10,000 samples, use two-stage approach
# Stage 1: Phase common variants (MAF > 0.1%)
shapeit5_phase_common \
    --input input.bcf \
    --map genetic_map.txt \
    --output phased_common.bcf \
    --thread 16 \
    --filter-maf 0.001

# Stage 2: Phase rare variants using common as scaffold
shapeit5_phase_rare \
    --input input.bcf \
    --scaffold phased_common.bcf \
    --map genetic_map.txt \
    --output phased_all.bcf \
    --thread 16
```

## Using Reference Panel for Better Phasing

```bash
# Beagle with 1000 Genomes reference
java -Xmx16g -jar beagle.jar \
    gt=study.vcf.gz \
    ref=1000GP.chr22.vcf.gz \
    map=plink.chr22.GRCh38.map \
    out=phased_with_ref \
    nthreads=8
```

## Input QC Before Phasing

```bash
# 1. Check for strand issues
bcftools +fixref input.vcf.gz -Oz -o fixed.vcf.gz -- -f reference.fa

# 2. Filter to biallelic SNPs
bcftools view -m2 -M2 -v snps fixed.vcf.gz -Oz -o biallelic.vcf.gz

# 3. Remove problematic regions
bcftools view -T ^exclude_regions.bed biallelic.vcf.gz -Oz -o filtered.vcf.gz

# 4. Check sample missingness
bcftools stats filtered.vcf.gz | grep "PSC"
```

## Output Format

Phased genotypes use `|` separator:
```
# Unphased
GT: 0/1  (unknown which allele on which haplotype)

# Phased
GT: 0|1  (REF on haplotype 1, ALT on haplotype 2)
GT: 1|0  (ALT on haplotype 1, REF on haplotype 2)
```
