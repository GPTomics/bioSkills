# Neoantigen Pipeline - Usage Guide

## Overview

This workflow identifies tumor-specific neoantigens from somatic mutations for personalized cancer vaccine design. It integrates HLA typing, MHC binding prediction, and multi-factor immunogenicity scoring to rank vaccine candidates.

## Prerequisites

```bash
pip install pvactools mhcflurry vatools pandas numpy matplotlib seaborn

mhcflurry-downloads fetch

conda install -c bioconda vep arcashla optitype samtools
```

**Required databases:**
- VEP cache for annotation
- IEDB tools (optional, for additional algorithms)

## Quick Start

Tell your AI agent what you want to do:
- "Find neoantigens from my somatic VCF for vaccine design"
- "Predict MHC binding for tumor mutations"
- "Rank neoantigen candidates by immunogenicity"
- "Run the neoantigen pipeline with my HLA types"

## Example Prompts

### Complete pipeline
> "Run neoantigen discovery on my tumor VCF with HLA-A*02:01,HLA-B*07:02"

> "Find vaccine candidates from my annotated somatic mutations"

### HLA typing
> "Determine HLA types from my tumor RNA-seq BAM"

> "Extract HLA alleles for neoantigen prediction"

### Binding prediction
> "Predict MHC Class I binding for my mutant peptides"

> "Find strong binders (<500nM) in my neoantigen candidates"

### Ranking
> "Score neoantigens by immunogenicity and expression"

> "Rank my neoantigen candidates for vaccine prioritization"

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| Somatic VCF | VCF (VEP-annotated) | Tumor somatic mutations |
| HLA types | String | Comma-separated 4-digit HLA alleles |
| Expression (optional) | TSV | Gene-level TPM from tumor RNA-seq |
| Tumor BAM (optional) | BAM | For HLA typing if types unknown |

## What the Agent Will Do

1. **HLA Typing** - Determines patient HLA alleles from RNA-seq (if not provided)
2. **VCF Annotation** - Adds protein consequences with VEP
3. **Binding Prediction** - Predicts peptide-MHC binding with multiple algorithms
4. **Neoantigen Calling** - Identifies tumor-specific peptides with pVACseq
5. **Immunogenicity Scoring** - Ranks candidates by binding, expression, VAF, and specificity
6. **Visualization** - Generates summary plots of candidate distribution

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| IC50 threshold | 500 nM | Strong binder cutoff |
| Epitope lengths | 8,9,10,11 | MHC-I peptide lengths |
| VAF minimum | 0.1 | Variant allele frequency filter |
| Expression minimum | 1 TPM | Gene expression filter |
| DAI threshold | 500 | Differential agretopicity for specificity |

## Tips

- **Expression data improves ranking**: Include tumor RNA-seq TPM when available
- **Use multiple algorithms**: MHCflurry + NetMHCpan gives more robust predictions
- **Consider Class II**: CD4+ T cell help improves vaccine efficacy
- **Clonal mutations first**: Prioritize high-VAF variants for broader tumor coverage; clonal neoantigen burden (McGranahan 2016) predicts ICI response better than total
- **Validate HLA typing**: Clinical-grade HLA typing is more reliable than computational; T1K (Song 2023 *Genome Res*) is the 2024-2026 all-rounder for class I + II + KIR from short read
- **Check HLA-LOH**: LOHHLA (Marty 2017) or DASH (Montesion 2021); somatic loss of HLA abolishes neoantigen presentation in ~17% of pan-cancer (>30% HNSCC / NSCLC / cervical); flag and exclude lost-allele predictions
- **Pair with TMB / MSI**: TMB-H pan-tumor (Vega 2021 panel-calibrated) and MSI-H / dMMR (KEYNOTE-016/164/158) are co-biomarkers; MSI-H supersedes TMB-H per Sha 2020

## Related Skills

- immunoinformatics/mhc-binding-prediction - MHCflurry parameters
- immunoinformatics/neoantigen-prediction - pVACtools details
- immunoinformatics/immunogenicity-scoring - Ranking algorithms
- clinical-databases/hla-typing - HLA typing + HLA-LOH detection
- clinical-databases/tumor-mutational-burden - TMB-H ICI biomarker
- clinical-databases/msi-detection - MSI-H / dMMR ICI biomarker
- clinical-databases/somatic-signatures - Mutational mechanism context
