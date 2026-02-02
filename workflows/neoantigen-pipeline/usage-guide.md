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
- **Clonal mutations first**: Prioritize high-VAF variants for broader tumor coverage
- **Validate HLA typing**: Clinical-grade HLA typing is more reliable than computational
