# Genotype Imputation Usage Guide

## Overview

Genotype imputation predicts untyped variants using haplotype patterns from a reference panel. This enables:
- Increased variant density for GWAS
- Fine-mapping of association signals
- Meta-analysis across different genotyping arrays
- Polygenic risk score calculation

## Imputation Workflow

```
1. QC and filter study data
2. Align to reference (strand, allele coding)
3. Phase haplotypes
4. Impute from reference panel
5. QC imputed variants (INFO score filter)
6. Association analysis
```

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

echo "Done! Output: $OUT_DIR/imputed.all.vcf.gz"
```

## Choosing a Reference Panel

| Panel | Size | Populations | Use For |
|-------|------|-------------|---------|
| 1000 Genomes | 2,504 | 26 global | General purpose |
| HRC | 32,470 | European-heavy | European studies |
| TOPMed | 97,256 | Diverse | Best accuracy |
| gnomAD | 76,156 | Diverse | Rare variants |

Match reference ancestry to your study population.

## Using Michigan Imputation Server

The easiest approach for most users:

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

## Interpreting Imputation Quality

### INFO/R2 Score
- Measures correlation between true and imputed genotypes
- Higher = better imputation
- Based on comparing expected vs observed variance

```bash
# Distribution of INFO scores
bcftools query -f '%INFO/DR2\n' imputed.vcf.gz | \
    awk '{if($1<0.3) low++; else if($1<0.8) med++; else high++}
    END {print "Low (<0.3):", low; print "Medium (0.3-0.8):", med; print "High (>0.8):", high}'
```

### Quality by MAF
Rare variants impute less accurately:

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load INFO scores and MAF
data = pd.read_csv('info_maf.txt', sep='\t', names=['CHR', 'POS', 'INFO', 'MAF'])

# Plot
plt.scatter(data['MAF'], data['INFO'], alpha=0.1, s=1)
plt.xlabel('Minor Allele Frequency')
plt.ylabel('INFO Score (R2)')
plt.title('Imputation Quality by MAF')
plt.savefig('info_by_maf.png')
```

## Using Imputed Data in GWAS

### With Dosages (Recommended)
```bash
# PLINK2
plink2 --vcf imputed.vcf.gz dosage=DS \
    --glm \
    --pheno phenotypes.txt \
    --covar covariates.txt \
    --out gwas_dosage

# REGENIE
regenie \
    --step 2 \
    --bgen imputed.bgen \
    --sample imputed.sample \
    --phenoFile phenotypes.txt \
    --covarFile covariates.txt \
    --out regenie_results
```

### With Hard Calls
```bash
# Convert dosage to hard call (loses uncertainty)
bcftools +dosage2gt imputed.vcf.gz -Oz -o imputed_gt.vcf.gz
```

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

## Troubleshooting

### Low INFO Scores
- Check reference panel ancestry match
- Verify strand alignment
- Remove low-quality genotyped variants

### High Missingness After Imputation
- Reference panel may not cover region
- Study variants not in reference

### Memory Errors
```bash
# Reduce chunk size
java -Xmx64g -jar beagle.jar window=20 ...

# Or process smaller regions
bcftools view -r chr22:1-25000000 input.vcf.gz ...
```
