# phasing-imputation

## Overview

Phase haplotypes and impute missing genotypes using reference panels. Essential for GWAS, population genetics, and integrating array and sequencing data. Covers Beagle, SHAPEIT, and reference panel management.

**Tool type:** cli | **Primary tools:** Beagle, SHAPEIT5, bcftools

## Skills

| Skill | Description |
|-------|-------------|
| haplotype-phasing | Phase genotypes into haplotypes with Beagle/SHAPEIT |
| genotype-imputation | Impute missing genotypes using reference panels |
| reference-panels | Work with 1000 Genomes, HRC, TOPMed reference panels |
| imputation-qc | Quality control of imputation results, INFO scores |

## Example Prompts

- "Phase my VCF file with Beagle"
- "Impute missing genotypes using 1000 Genomes"
- "Download and prepare reference panel for imputation"
- "Filter imputed variants by INFO score"
- "Run SHAPEIT5 on my GWAS data"
- "Check imputation quality metrics"
- "Convert reference panel to Beagle format"
- "Phase and impute chromosome by chromosome"
- "Compare phasing accuracy"
- "Prepare data for imputation server"

## Requirements

```bash
# Beagle
wget https://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar
# Run with: java -jar beagle.jar ...

# SHAPEIT5 (for large biobank data)
conda install -c bioconda shapeit5

# bcftools (for VCF manipulation)
conda install -c bioconda bcftools

# Reference panels (download separately)
# 1000 Genomes, HRC, TOPMed
```

## Related Skills

- **variant-calling** - Generate input VCF files, VCF processing
- **population-genetics** - Downstream analysis (PCA, GWAS)
- **genome-intervals** - Region-based filtering
