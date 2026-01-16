# PLINK Basics - Usage Guide

## Overview

PLINK is the standard tool for population genetic analysis. PLINK 2.0 is faster and more memory-efficient but doesn't support all legacy formats. Use PLINK 1.9 for format conversion from PED/MAP, then PLINK 2.0 for analysis.

## When to Use This Skill

- Converting VCF to PLINK format for analysis
- Quality control filtering before GWAS
- Merging datasets from different sources
- Extracting specific samples or variants
- Basic allele frequency calculations

## Installation

```bash
# PLINK 1.9
conda install -c bioconda plink

# PLINK 2.0
conda install -c bioconda plink2
```

## Standard QC Workflow

### 1. Convert VCF to PLINK

```bash
plink2 --vcf data.vcf.gz --double-id --make-bed --out data
```

### 2. Check Initial Statistics

```bash
# Sample and variant counts
wc -l data.fam data.bim

# Missing rates
plink2 --bfile data --missing --out data_missing
```

### 3. Apply QC Filters

```bash
plink2 --bfile data \
    --maf 0.01 \
    --geno 0.05 \
    --mind 0.05 \
    --hwe 1e-6 \
    --make-bed --out data_qc
```

### 4. Report Filtering

```bash
echo "Before QC:"
wc -l data.fam data.bim

echo "After QC:"
wc -l data_qc.fam data_qc.bim
```

## Choosing Thresholds

| Filter | Conservative | Standard | Lenient |
|--------|--------------|----------|---------|
| MAF | 0.05 | 0.01 | 0.001 |
| Geno | 0.02 | 0.05 | 0.10 |
| Mind | 0.02 | 0.05 | 0.10 |
| HWE | 1e-4 | 1e-6 | 1e-10 |

## Common Issues

### ID Mismatch

VCF sample IDs may need parsing:
```bash
# Double ID (same FID and IID)
plink2 --vcf input.vcf.gz --double-id --make-bed --out output
```

### Allele Code Issues

```bash
# Handle non-ACGT alleles
plink2 --bfile input --snps-only just-acgt --make-bed --out output
```

### Duplicate IDs

```bash
# Check for duplicates
awk '{print $2}' data.bim | sort | uniq -d

# Remove duplicates
plink2 --bfile input --rm-dup force-first --make-bed --out output
```

## Resources

- [PLINK 1.9 Documentation](https://www.cog-genomics.org/plink/1.9/)
- [PLINK 2.0 Documentation](https://www.cog-genomics.org/plink/2.0/)
