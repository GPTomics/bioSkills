# causal-genomics

## Overview

Infer causal relationships from genetic association data using Mendelian randomization, colocalization, and mediation analysis.

**Tool type:** r | **Primary tools:** TwoSampleMR, coloc, susieR, mediation, MR-PRESSO

## Skills

| Skill | Description |
|-------|-------------|
| mendelian-randomization | Causal inference from GWAS with TwoSampleMR |
| colocalization-analysis | Bayesian colocalization with coloc |
| mediation-analysis | Causal mediation with the mediation package |
| fine-mapping | Credible set construction with susieR, FINEMAP |
| pleiotropy-detection | Pleiotropy diagnostics with MR-PRESSO, MR-Egger |

## Example Prompts

- "Test whether BMI causally affects type 2 diabetes using GWAS summary statistics"
- "Check if this GWAS signal and eQTL share a causal variant"
- "Does gene expression mediate the effect of this SNP on disease risk?"
- "Fine-map this GWAS locus to a 95% credible set"
- "Run sensitivity analyses for my MR results to check for pleiotropy"

## Requirements

```bash
# R packages from CRAN
install.packages(c('remotes', 'mediation'))

# TwoSampleMR (from GitHub)
remotes::install_github('MRCIEU/TwoSampleMR')

# coloc
install.packages('coloc')

# susieR
install.packages('susieR')

# MR-PRESSO
remotes::install_github('rondolab/MR-PRESSO')
```

## Related Skills

- **population-genetics** - GWAS, LD, and association testing
- **clinical-databases** - ClinVar, gnomAD variant annotation
- **pathway-analysis** - GO enrichment and GSEA for causal gene sets
- **differential-expression** - eQTL data for colocalization and mediation
