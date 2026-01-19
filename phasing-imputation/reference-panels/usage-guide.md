# Reference Panels Usage Guide

## Overview

Reference panels provide haplotype information from well-characterized samples. They are essential for:
- Haplotype phasing (using LD patterns)
- Genotype imputation (predicting untyped variants)
- Quality control (strand and allele alignment)

## Choosing a Reference Panel

### By Ancestry
Match the reference panel to your study population:

| Study Population | Recommended Panel |
|-----------------|-------------------|
| European | HRC, 1000G EUR, UK10K |
| African | 1000G AFR, TOPMed |
| East Asian | 1000G EAS |
| South Asian | 1000G SAS |
| Mixed/diverse | TOPMed, 1000G ALL |
| Latino/Admixed | TOPMed, 1000G AMR |

### By Application
| Application | Recommended |
|-------------|-------------|
| General GWAS | 1000G or HRC |
| Rare variants | TOPMed |
| Fine-mapping | TOPMed or HRC |
| Non-European | 1000G superpopulation |

## Setting Up 1000 Genomes

```bash
#!/bin/bash
# Download and prepare 1000 Genomes Phase 3 (GRCh38)

BASE_DIR=reference_panels/1000GP
mkdir -p $BASE_DIR

# Download high-coverage data
FTP="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased"

for chr in {1..22}; do
    echo "Downloading chromosome $chr..."
    wget -P $BASE_DIR ${FTP}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz
    wget -P $BASE_DIR ${FTP}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz.tbi
done

# Download sample information
wget -P $BASE_DIR http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt

# Create population subsets
cd $BASE_DIR
awk '$7=="EUR" {print $2}' 20130606_g1k_3202_samples_ped_population.txt > EUR_samples.txt
awk '$7=="AFR" {print $2}' 20130606_g1k_3202_samples_ped_population.txt > AFR_samples.txt
awk '$7=="EAS" {print $2}' 20130606_g1k_3202_samples_ped_population.txt > EAS_samples.txt
```

## Preparing Reference for Imputation

```bash
# Standard preparation workflow

# 1. Filter to biallelic SNPs
bcftools view -m2 -M2 -v snps reference.vcf.gz -Oz -o ref_biallelic.vcf.gz

# 2. Remove rare variants if needed (optional, saves memory)
bcftools view -q 0.001:minor ref_biallelic.vcf.gz -Oz -o ref_maf001.vcf.gz

# 3. Set consistent IDs
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' ref_maf001.vcf.gz -Oz -o ref_final.vcf.gz

# 4. Index
bcftools index ref_final.vcf.gz
```

## Downloading Genetic Maps

```bash
# Beagle format (GRCh38)
mkdir -p genetic_maps
cd genetic_maps

wget https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip
unzip plink.GRCh38.map.zip

# Files: plink.chr1.GRCh38.map, plink.chr2.GRCh38.map, etc.
# Format: position COMBINED_rate(cM/Mb) Genetic_Map(cM)
```

## Aligning Study Data to Reference

Before imputation, ensure study data matches reference:

```bash
# 1. Check reference genome build
bcftools view -h study.vcf.gz | grep "##reference"
bcftools view -h reference.vcf.gz | grep "##reference"

# 2. Fix strand and allele issues
bcftools +fixref study.vcf.gz -Oz -o study_fixed.vcf.gz -- \
    -f genome.fa \
    -i reference.vcf.gz \
    -m flip

# 3. Extract overlapping sites
bcftools isec -n=2 -w1 \
    study_fixed.vcf.gz \
    reference.vcf.gz \
    -Oz -o study_for_imputation.vcf.gz
```

## Liftover Between Builds

```bash
# GRCh37 to GRCh38 using Picard
# Download chain file from UCSC

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

java -jar picard.jar LiftoverVcf \
    I=study_hg19.vcf.gz \
    O=study_hg38.vcf.gz \
    CHAIN=hg19ToHg38.over.chain.gz \
    R=hg38.fa \
    REJECT=rejected.vcf

# Check rejection rate
wc -l rejected.vcf  # Should be <5% typically
```

## Using Imputation Servers

For convenience, use web-based imputation servers:

### Michigan Imputation Server
```bash
# Prepare files
for chr in {1..22}; do
    bcftools view -r chr${chr} study.vcf.gz -Oz -o upload/study.chr${chr}.vcf.gz
done

# Upload to https://imputationserver.sph.umich.edu
# Select:
# - Reference Panel: HRC r1.1 2016 or TOPMed r2
# - Population: matching your study
# - Phasing: Eagle (default)
```

### TOPMed Imputation Server
```bash
# Same preparation, upload to:
# https://imputation.biodatacatalyst.nhlbi.nih.gov

# TOPMed has best coverage for rare variants
# Requires dbGaP authorization for some uses
```

## Quality Check Reference Panel

```python
import subprocess
import pandas as pd

def check_reference_panel(vcf_path):
    '''Basic QC of reference panel.'''
    # Sample count
    samples = subprocess.check_output(f'bcftools query -l {vcf_path} | wc -l', shell=True)
    print(f'Samples: {samples.decode().strip()}')

    # Variant count
    variants = subprocess.check_output(f'bcftools view -H {vcf_path} | wc -l', shell=True)
    print(f'Variants: {variants.decode().strip()}')

    # Chromosomes
    chroms = subprocess.check_output(f'bcftools index -s {vcf_path}', shell=True)
    print(f'Chromosomes:\n{chroms.decode()}')

    # Check phasing
    gt_sample = subprocess.check_output(
        f'bcftools query -f "[%GT\\n]" {vcf_path} | head -100', shell=True)
    phased = gt_sample.decode().count('|')
    total = gt_sample.decode().count('\n')
    print(f'Phased genotypes: {phased}/{total} ({100*phased/total:.1f}%)')

check_reference_panel('reference.vcf.gz')
```

## Storage Requirements

| Panel | Compressed Size | Uncompressed |
|-------|-----------------|--------------|
| 1000G Phase 3 | ~50 GB | ~500 GB |
| HRC | ~100 GB | ~1 TB |
| TOPMed | ~500 GB | ~5 TB |

Plan storage accordingly, especially for per-chromosome files.
