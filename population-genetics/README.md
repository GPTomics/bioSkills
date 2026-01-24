# population-genetics

## Overview

Population genetic analysis using PLINK, Admixture, and scikit-allel. Covers GWAS, population structure analysis, selection statistics, and linkage disequilibrium calculations.

**Tool type:** mixed | **Primary tools:** PLINK 1.9/2.0, FlashPCA2, Admixture, scikit-allel

## Skills

| Skill | Description |
|-------|-------------|
| plink-basics | File formats, conversion, QC filtering |
| association-testing | GWAS, case-control, quantitative trait association |
| population-structure | PCA (PLINK2/FlashPCA2), admixture analysis, MDS |
| scikit-allel-analysis | Python population genetics |
| selection-statistics | Fst, Tajima's D, iHS, XP-EHH |
| linkage-disequilibrium | LD calculations, pruning, haplotype blocks |

## Example Prompts

- "Convert VCF to PLINK binary format"
- "Filter variants by MAF and missingness"
- "Apply QC filters to my GWAS data"
- "Run GWAS on my case-control data"
- "Test association for a quantitative trait"
- "Create a Manhattan plot of my GWAS results"
- "Perform PCA for population structure"
- "Run fast PCA on biobank-scale data with FlashPCA2"
- "Run admixture analysis with K=3"
- "Plot admixture results as a bar chart"
- "Calculate Fst between populations"
- "Compute Tajima's D in windows"
- "Calculate iHS for selection scan"
- "Prune SNPs by linkage disequilibrium"
- "Calculate LD between two SNPs"
- "Identify haplotype blocks"

## Requirements

```bash
# PLINK 1.9 and 2.0
conda install -c bioconda plink plink2

# Admixture
conda install -c bioconda admixture

# FlashPCA2 (for large datasets)
conda install -c bioconda flashpca

# scikit-allel
pip install scikit-allel
```

## Related Skills

- **variant-calling** - VCF generation, filtering, and manipulation
- **phasing-imputation** - Haplotype phasing for population studies
- **alignment-files** - BAM processing before variant calling
