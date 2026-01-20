# Genotype Imputation - Usage Guide

## Overview
Genotype imputation predicts untyped variants using haplotype patterns from a reference panel, enabling increased variant density for GWAS, fine-mapping, and polygenic risk score calculation.

## Prerequisites
```bash
# Beagle (download jar)
wget https://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar

# bcftools for data manipulation
conda install -c bioconda bcftools

# PLINK2 for downstream analysis
conda install -c bioconda plink2

# Genetic maps
wget https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip
unzip plink.GRCh38.map.zip -d genetic_maps/

# Reference panel (e.g., 1000 Genomes)
# Download from http://ftp.1000genomes.ebi.ac.uk/
```

## Quick Start
Tell your AI agent what you want to do:
- "Impute missing genotypes in my VCF using 1000 Genomes reference"
- "Set up an imputation pipeline for all chromosomes"
- "Prepare my data for the Michigan Imputation Server"

## Example Prompts
### Basic Imputation
> "Impute genotypes in study.vcf.gz using 1000 Genomes Phase 3 as reference panel"

### Full Pipeline
> "Run a complete imputation pipeline: align to reference, phase, impute, and filter by INFO score"

### Server Preparation
> "Prepare my data for upload to the Michigan Imputation Server"

### Quality Filtering
> "Filter my imputed VCF to keep only variants with INFO score > 0.8 and MAF > 0.01"

## What the Agent Will Do
1. Align study data to reference panel (strand and allele coding)
2. Phase haplotypes if not already phased
3. Impute genotypes from reference panel per chromosome
4. Filter imputed variants by INFO score
5. Merge chromosomes into final output

## Tips
- Match reference panel ancestry to your study population for best results
- TOPMed provides best accuracy for rare variants
- Michigan Imputation Server is the easiest approach for most users
- Always filter by INFO score (typically > 0.3 for GWAS)
- Use dosages (not hard calls) in association analysis
- Rare variants impute less accurately - check INFO by MAF

## Choosing a Reference Panel

| Panel | Size | Populations | Use For |
|-------|------|-------------|---------|
| 1000 Genomes | 2,504 | 26 global | General purpose |
| HRC | 32,470 | European-heavy | European studies |
| TOPMed | 97,256 | Diverse | Best accuracy |
| gnomAD | 76,156 | Diverse | Rare variants |

## Complete Pipeline

```bash
#!/bin/bash
set -euo pipefail

STUDY=study.vcf.gz
REF_DIR=reference_panels/1000GP_phase3
OUT_DIR=imputed
THREADS=8

mkdir -p $OUT_DIR

for chr in {1..22}; do
    echo "Processing chromosome $chr..."

    # Extract chromosome
    bcftools view -r chr${chr} $STUDY -Oz -o $OUT_DIR/study.chr${chr}.vcf.gz

    # Align to reference
    bcftools +fixref $OUT_DIR/study.chr${chr}.vcf.gz \
        -Oz -o $OUT_DIR/fixed.chr${chr}.vcf.gz -- \
        -f reference.fa -m flip

    # Phase and impute together (Beagle does both)
    java -Xmx32g -jar beagle.jar \
        gt=$OUT_DIR/fixed.chr${chr}.vcf.gz \
        ref=${REF_DIR}/1000GP.chr${chr}.vcf.gz \
        map=genetic_maps/plink.chr${chr}.GRCh38.map \
        out=$OUT_DIR/imputed.chr${chr} \
        gp=true \
        nthreads=$THREADS

    # Filter by imputation quality
    bcftools view -i 'INFO/DR2 > 0.3' \
        $OUT_DIR/imputed.chr${chr}.vcf.gz \
        -Oz -o $OUT_DIR/imputed.chr${chr}.filtered.vcf.gz

    bcftools index $OUT_DIR/imputed.chr${chr}.filtered.vcf.gz
done

# Merge chromosomes
bcftools concat $OUT_DIR/imputed.chr*.filtered.vcf.gz \
    -Oz -o $OUT_DIR/imputed.all.vcf.gz
bcftools index $OUT_DIR/imputed.all.vcf.gz
```

## Using Michigan Imputation Server

```bash
# 1. Prepare files (VCF per chromosome)
for chr in {1..22}; do
    bcftools view -r chr${chr} study.vcf.gz -Oz -o study.chr${chr}.vcf.gz
done

# 2. Upload to imputationserver.sph.umich.edu
#    - Select reference panel
#    - Choose population
#    - Enable phasing

# 3. Download results
#    - Imputed VCF with dosages
#    - INFO score file
#    - QC report
```

## Interpreting INFO Scores

| R2 Score | Quality | Recommendation |
|----------|---------|----------------|
| > 0.9 | Excellent | All analyses |
| 0.8 - 0.9 | Good | Fine-mapping, PRS |
| 0.5 - 0.8 | Moderate | GWAS |
| 0.3 - 0.5 | Low | GWAS discovery only |
| < 0.3 | Poor | Exclude |

## Post-Imputation QC

```bash
# 1. Filter by INFO
bcftools view -i 'INFO/DR2 > 0.3' imputed.vcf.gz -Oz -o filtered.vcf.gz

# 2. Filter by MAF
bcftools view -i 'MAF > 0.01' filtered.vcf.gz -Oz -o common.vcf.gz

# 3. Remove palindromic SNPs (optional)
bcftools view -e 'REF="A" && ALT="T" || REF="T" && ALT="A" || REF="G" && ALT="C" || REF="C" && ALT="G"' \
    common.vcf.gz -Oz -o no_palindrome.vcf.gz

# 4. Hardy-Weinberg filter
plink2 --vcf common.vcf.gz --hwe 1e-6 --make-pgen --out qc_passed
```

## Using Imputed Data in GWAS

```bash
# With dosages (recommended)
plink2 --vcf imputed.vcf.gz dosage=DS \
    --glm \
    --pheno phenotypes.txt \
    --covar covariates.txt \
    --out gwas_dosage

# Convert to hard calls if needed (loses uncertainty)
bcftools +dosage2gt imputed.vcf.gz -Oz -o imputed_gt.vcf.gz
```
