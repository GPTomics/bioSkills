# phasing-imputation

## Overview

Statistical haplotype phasing and reference-panel genotype imputation for array and sequence data, with Beagle, SHAPEIT5, IMPUTE5, Minimac4, and GLIMPSE2. Decision-grade framing: phasing and imputation are one Li-Stephens copying HMM, so phase is a statistical inference (not a measurement, and distinct from read-backed phasing) and imputation's honest output is a dosage with a self-estimated quality (R2/DR2/INFO) that must be filtered before GWAS; the panel is the prior, so its ancestry must match the target; and strand and genome-build alignment to the panel is mandatory or accuracy silently collapses. Read-backed single-sample phasing lives in long-read-sequencing.

**Tool type:** cli | **Primary tools:** Beagle, SHAPEIT5, bcftools

## Skills

| Skill | Description |
|-------|-------------|
| reference-panels | Select and prepare 1000G/HRC/TOPMed panels; ancestry-match, build, strand, server access |
| haplotype-phasing | Statistical phasing with SHAPEIT5/Eagle2/Beagle; switch error; the substrate for imputation |
| genotype-imputation | Impute untyped variants with Beagle/Minimac4/IMPUTE5/GLIMPSE2; dosages |
| imputation-qc | Filter imputed variants by DR2/R2/INFO and MAF; masked-truth accuracy |

## Method Selection

| Scenario | Method | Skill |
|----------|--------|-------|
| Common-variant GWAS, well-paneled ancestry | array + impute against a matched panel | genotype-imputation |
| Rare variants / under-represented ancestry | low-coverage WGS + GLIMPSE2 | genotype-imputation |
| Biobank-scale phasing (100k+), rare variants | SHAPEIT5 | haplotype-phasing |
| Moderate cohort, one tool for phase + impute | Beagle | haplotype-phasing, genotype-imputation |
| Diverse / admixed target | TOPMed panel (server) | reference-panels |
| Predominantly European | HRC (server) or 1000G | reference-panels |
| Diverse and downloadable (local) | HGDP+1kGP (gnomAD) | reference-panels |
| Single individual, long reads | (read-backed) | long-read-sequencing/haplotype-phasing |

## Example Prompts

- "Phase my VCF file with Beagle"
- "Impute missing genotypes against TOPMed"
- "Should I use an array or low-coverage WGS for my admixed cohort?"
- "Select and prepare a reference panel for my target ancestry"
- "Filter imputed variants by R2 and MAF before GWAS"
- "Run SHAPEIT5 rare-variant phasing on my biobank data"
- "Check imputation quality stratified by allele frequency"
- "Align strand and genome build to the panel before imputing"
- "Phase and impute chromosome by chromosome"
- "Validate imputation accuracy by masking typed genotypes"

## Requirements

```bash
# Beagle (Java jar)
wget https://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar
# Run with: java -jar beagle.jar ...

# SHAPEIT5 (biobank phasing), Minimac4/IMPUTE5/GLIMPSE2 (imputation), bcftools
conda install -c bioconda shapeit5 minimac4 bcftools

# Reference panels (download separately; HRC/TOPMed are server-only)
# 1000 Genomes, HGDP+1kGP (downloadable); HRC, TOPMed (imputation servers)
```

## Related Skills

- **variant-calling** - Generate and normalize input VCF files
- **population-genetics** - Association testing, population structure, LD
- **long-read-sequencing** - Read-backed / molecular single-sample phasing
- **workflows** - gwas-pipeline orchestrates QC -> phase -> impute -> associate
