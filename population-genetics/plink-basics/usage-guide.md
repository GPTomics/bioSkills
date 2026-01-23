# PLINK Basics - Usage Guide

## Overview

PLINK is the standard tool for population genetic analysis, handling format conversion, quality control, and basic statistics. PLINK 2.0 is faster and more memory-efficient but doesn't support all legacy formats.

## Prerequisites

```bash
# PLINK 1.9
conda install -c bioconda plink

# PLINK 2.0
conda install -c bioconda plink2
```

## Quick Start

Tell your AI agent what you want to do:
- "Convert my VCF file to PLINK format"
- "Run quality control on my GWAS data"
- "Filter out low-quality variants and samples"
- "Merge these two PLINK datasets"
- "Extract samples from a specific population"

## Example Prompts

### Format Conversion
> "Convert data.vcf.gz to PLINK binary format"

> "Convert my PED/MAP files to binary BED format"

> "Export my PLINK data back to VCF"

### Quality Control
> "Run standard QC filtering on my GWAS data with MAF 0.01, genotype missingness 5%, and sample missingness 5%"

> "Apply strict QC filters for population structure analysis"

> "Check for samples with high missingness rates"

### Data Management
> "Merge these three PLINK datasets into one"

> "Extract only the European samples from my data"

> "Remove all variants not in my keep list"

> "Calculate allele frequencies for each population"

## What the Agent Will Do

1. Assess input data format (VCF, PED/MAP, or binary)
2. Determine appropriate PLINK version (1.9 for legacy formats, 2.0 for analysis)
3. Run quality control or conversion commands
4. Report statistics before and after filtering
5. Verify output file integrity

## Tips

- Use `--double-id` when converting VCF to handle sample ID parsing
- PLINK 2.0's `--pfile` is faster than `--bfile` for large datasets
- Always check variant/sample counts before and after QC
- Use `--snps-only just-acgt` to remove indels and non-standard alleles
- For merging datasets, ensure consistent chromosome naming and strand orientation

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
