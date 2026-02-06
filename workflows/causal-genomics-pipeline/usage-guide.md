# Causal Genomics Pipeline - Usage Guide

## Overview
Complete post-GWAS causal inference workflow from summary statistics through Mendelian randomization, sensitivity analysis, colocalization, fine-mapping, and mediation analysis. Integrates multiple complementary approaches to triangulate causal evidence between exposures and outcomes.

## Prerequisites
```r
# Core MR packages
install.packages('remotes')
remotes::install_github('MRCIEU/TwoSampleMR')
remotes::install_github('rondolab/MR-PRESSO')

# Colocalization
install.packages('coloc')

# Fine-mapping
install.packages('susieR')

# Multivariable MR
install.packages('MendelianRandomization')
```

**Input data:**
- GWAS summary statistics for exposure and outcome (TSV with SNP, BETA, SE, A1, A2, EAF, P columns)
- LD matrix for fine-mapping loci (optional, from reference panel)
- For colocalization: full summary statistics for overlapping loci

## Quick Start
Tell your AI agent what you want to do:
- "Run Mendelian randomization from my GWAS summary statistics"
- "Perform post-GWAS causal inference with MR and colocalization"
- "Test whether my exposure is causally related to the outcome"
- "Fine-map the causal variant at my GWAS locus"

## Example Prompts

### Full Pipeline
> "I have GWAS summary statistics for LDL cholesterol and coronary heart disease. Run the complete causal inference pipeline."

> "Perform MR, sensitivity analysis, and colocalization to test if BMI causally affects type 2 diabetes."

### Specific Steps
> "Select strong instruments from my exposure GWAS and check F-statistics."

> "Run MR-PRESSO to detect outlier instruments and correct the MR estimate."

> "Colocalize my GWAS locus with an eQTL dataset to find the shared causal variant."

> "Fine-map my GWAS locus with SuSiE and report 95% credible sets."

### Mediation
> "Test if CRP mediates the causal effect of BMI on coronary disease using MVMR."

> "Run network MR to identify the causal pathway between my exposures and outcome."

## What the Agent Will Do
1. Load GWAS summary statistics for exposure and outcome
2. Select independent instruments (LD clumping, F-stat > 10)
3. Run MR across multiple methods (IVW, Egger, weighted median, weighted mode)
4. Perform sensitivity analysis (MR-PRESSO, Egger intercept, Steiger test)
5. Run colocalization at top loci (coloc.abf)
6. Fine-map causal variants with SuSiE (95% credible sets)
7. Optionally test mediation with multivariable MR
8. Summarize triangulated causal evidence

## Tips
- Always use multiple MR methods; if they disagree, investigate pleiotropy
- F-statistic > 10 is minimum; > 30 preferred to avoid weak instrument bias
- MR-PRESSO global test p > 0.05 means no significant horizontal pleiotropy
- PP.H4 > 0.8 is strong colocalization evidence; 0.5-0.8 is suggestive
- Colocalization priors (p12) matter: 1e-5 is standard, 5e-6 is conservative
- Fine-mapping requires an accurate LD matrix matching the study population
- Always check Steiger directionality to rule out reverse causation
- For mediation, ensure the mediator GWAS is adequately powered

## Related Skills
- causal-genomics/mendelian-randomization - MR method details
- causal-genomics/colocalization-analysis - Coloc implementation
- causal-genomics/fine-mapping - SuSiE and FINEMAP details
- causal-genomics/pleiotropy-detection - MR-PRESSO and sensitivity tests
- causal-genomics/mediation-analysis - MVMR and network MR
- population-genetics/association-testing - Upstream GWAS methods
