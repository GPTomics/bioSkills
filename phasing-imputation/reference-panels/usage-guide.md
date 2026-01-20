# Reference Panels - Usage Guide

## Overview
Reference panels provide haplotype information from well-characterized samples, essential for haplotype phasing, genotype imputation, and quality control alignment.

## Prerequisites
```bash
# bcftools for panel manipulation
conda install -c bioconda bcftools

# Picard for liftover
conda install -c bioconda picard

# wget for downloading panels
# (usually pre-installed)

# Storage: 50GB-500GB depending on panel
```

## Quick Start
Tell your AI agent what you want to do:
- "Download and set up 1000 Genomes reference panel for imputation"
- "Create population-specific subsets from 1000 Genomes"
- "Liftover my data from GRCh37 to GRCh38"

## Example Prompts
### Panel Setup
> "Download 1000 Genomes Phase 3 high-coverage data and prepare it for imputation"

### Population Subset
> "Extract European samples from the 1000 Genomes reference panel"

### Data Alignment
> "Align my study VCF to the reference panel, fixing strand and allele issues"

### Genome Build Conversion
> "Liftover my GRCh37 VCF to GRCh38 to match the reference panel"

## What the Agent Will Do
1. Download reference panel files for all chromosomes
2. Extract population-specific samples if needed
3. Prepare reference for imputation (biallelic SNPs, consistent IDs)
4. Align study data to reference (strand correction)
5. Verify alignment and compatibility

## Tips
- Match reference panel ancestry to your study population
- TOPMed has best coverage for rare variants
- Always verify genome build matches between study and reference
- Plan storage: 1000G ~50GB, HRC ~100GB, TOPMed ~500GB
- Use imputation servers (Michigan, TOPMed) for convenience
- Keep population subset files for ancestry-matched imputation

## Choosing a Reference Panel

### By Ancestry
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

## Aligning Study Data to Reference

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
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

java -jar picard.jar LiftoverVcf \
    I=study_hg19.vcf.gz \
    O=study_hg38.vcf.gz \
    CHAIN=hg19ToHg38.over.chain.gz \
    R=hg38.fa \
    REJECT=rejected.vcf

# Check rejection rate (should be <5%)
wc -l rejected.vcf
```

## Using Imputation Servers

### Michigan Imputation Server
```bash
# Prepare files
for chr in {1..22}; do
    bcftools view -r chr${chr} study.vcf.gz -Oz -o upload/study.chr${chr}.vcf.gz
done

# Upload to https://imputationserver.sph.umich.edu
# Select: Reference Panel (HRC r1.1 or TOPMed r2), Population, Phasing (Eagle)
```

### TOPMed Imputation Server
```bash
# Same preparation, upload to:
# https://imputation.biodatacatalyst.nhlbi.nih.gov

# TOPMed has best coverage for rare variants
# May require dbGaP authorization
```

## Storage Requirements

| Panel | Compressed Size | Uncompressed |
|-------|-----------------|--------------|
| 1000G Phase 3 | ~50 GB | ~500 GB |
| HRC | ~100 GB | ~1 TB |
| TOPMed | ~500 GB | ~5 TB |
