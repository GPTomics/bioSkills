# Association Testing - Usage Guide

## Overview

Single-variant GWAS fits one regression per variant with an additive genotype dosage as the predictor, and the statistic is honest only insofar as genotype is independent of unmodeled phenotype drivers after the chosen covariates and random effects. The engine is set by sample structure and case:control imbalance, not taste: a PC-adjusted GLM (plink2) for unrelated samples whose structure top PCs capture, a linear mixed model (GEMMA, BOLT-LMM, SAIGE, regenie) the moment relatedness or fine-scale structure is present, and SPA or Firth when case:control imbalance and low minor-allele count break the score/Wald tail. Genomic inflation above 1 is expected under a polygenic trait and is mostly true signal; the LDSC intercept, not lambda, is the confounding diagnostic.

## Prerequisites

- plink2 installed: `conda install -c bioconda plink2`
- For mixed models: SAIGE, regenie, BOLT-LMM, or GEMMA (`conda install -c bioconda saige regenie gemma`)
- For plots and lambda: `pip install pandas matplotlib scipy numpy`
- Post-QC genotypes (see plink-basics) with phenotype coding verified, plus PCs (see population-structure)
- Conceptual prerequisites and big notes:
  - Genomic inflation (lambda) above 1 is EXPECTED under polygenicity and mostly true signal; do not divide statistics by lambda. The LDSC intercept separates confounding from polygenicity.
  - PC covariates absorb continuous ancestry but CANNOT remove relatedness (a covariance structure); a related or fine-scale-structured sample requires a linear mixed model.
  - Leave-one-chromosome-out (LOCO) is required in any LMM: a GRM that includes the candidate chromosome deflates power (proximal contamination).
  - SPA (SAIGE) and Firth (plink2 firth-fallback, regenie) exist because plain score/Wald tests are anti-conservative or collapse at extreme case:control imbalance and low MAC.
  - HWE filtering for QC is controls-only; filtering cases removes true non-additive associations.
  - Regress on imputed DOSAGES, not hard calls, and carry effect-allele/strand bookkeeping (CHR, POS, EA, OA, EAF) for any downstream meta-analysis or PRS.
  - Rare-variant gene-based aggregation (burden, SKAT, SKAT-O, ACAT) is a SEPARATE skill (rare-variant-association); this skill is single-variant only.

## Quick Start

Tell your AI agent what you want to do:
- "Run a single-variant GWAS for my case-control phenotype with PCs as covariates"
- "My sample has related individuals, pick the right association model"
- "My case:control ratio is 1:50, which test stays calibrated"
- "Regress on imputed dosages instead of hard calls"
- "Compute lambda and tell me whether it means confounding or polygenicity"
- "Run a biobank-scale GWAS with regenie step 1 and step 2"

## Example Prompts

### Single-variant GWAS
> "Run a logistic GWAS on my QC'd case-control data with the first 10 PCs, age, and sex as covariates, and write a Manhattan and QQ plot."

> "Test association with my quantitative trait on imputed dosages, reporting A1 frequency and imputation R2."

### Choosing the model
> "I have a cohort with cryptic relatedness and a binary trait with a 1:80 case:control ratio. Which association engine should I use and why?"

> "My quantitative-trait sample is 200,000 people with some families. Set up a mixed-model GWAS with LOCO."

### Calibration and confounding
> "My lambda is 1.18. Is that confounding or polygenicity, and what should I check?"

> "Set up SPA or Firth so my rare, imbalanced variants are not anti-conservative."

### Bookkeeping for downstream use
> "Prepare my GWAS sumstats with harmonized effect alleles so they can be meta-analyzed and used for a PRS."

## What the Agent Will Do

1. Confirm the genotypes are post-QC, the phenotype coding is correct, and PCs are available.
2. Diagnose sample structure and case:control balance to choose the engine: PC-adjusted GLM, an LMM, or an LMM with SPA/Firth.
3. Run association on dosages where imputed, with LOCO enabled for any mixed model.
4. Apply Firth or SPA when case:control imbalance or low MAC would break the score/Wald tail.
5. Compute lambda and (where sumstats allow) the LDSC intercept, and interpret inflation as polygenicity vs confounding rather than auto-correcting.
6. Report effect sizes with harmonized effect-allele/strand bookkeeping, flag winner's curse for downstream use, and route rare-variant aggregation to rare-variant-association.

## Tips

- Lambda above 1 is mostly polygenic signal under a real trait; do NOT genomic-control on lambda. Use the LDSC intercept (intercept ~ 1 with high lambda means polygenicity) and rescale to lambda_1000 before comparing studies of different N.
- The moment relatedness or fine-scale structure is present, switch to an LMM. Adding more PCs cannot remove a covariance structure.
- Keep LOCO on in every mixed model; a GRM that includes the candidate chromosome silently deflates power.
- Use Firth or SPA whenever case:control is more extreme than ~1:10 or MAC is low; a plain score/Wald test misses real rare-variant signal or fills the significant tail with artifacts.
- Compute HWE QC in controls only; filtering cases can remove a true non-additive disease variant.
- Do not adjust for a heritable covariate (BMI, smoking) without care; it can open a collider path and manufacture direction-flipped hits.
- Carry CHR, POS, effect allele, other allele, and effect-allele frequency in sumstats; resolve A/T and C/G palindromes by frequency or drop them.

## Related Skills

- plink-basics - QC, phenotype encoding, and the fileset that enters association
- population-structure - PCA covariates for stratification control
- linkage-disequilibrium - LD pruning before PCA and clumping of GWAS hits
- rare-variant-association - gene-based aggregation (burden, SKAT, SKAT-O, ACAT) below the single-variant MAC floor
- causal-genomics/fine-mapping - from an associated locus to a credible set of causal variants
- causal-genomics/mendelian-randomization - GWAS variants as instruments for causal inference
- clinical-databases/polygenic-risk - PRS built from association sumstats
- phasing-imputation/genotype-imputation - imputed dosages that enter --glm
