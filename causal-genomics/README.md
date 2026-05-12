# causal-genomics

## Overview

Infer causal relationships from genetic association data using Mendelian randomization, colocalization, fine-mapping, mediation, heritability partitioning, transcriptome-wide association, and drug-target validation with cis-pQTL MR.

**Tool type:** mixed | **Primary tools:** TwoSampleMR, MendelianRandomization, coloc, susieR, mediation, CMAverse, MR-PRESSO, CAUSE, LHC-MR, LDSC, LDAK, HDL, LAVA, FUSION, MetaXcan, FOCUS, MAGMA, PoPS, GenomicSEM, MTAG

## Skills

| Skill | Description |
|-------|-------------|
| mendelian-randomization | Causal inference from GWAS with TwoSampleMR, MendelianRandomization, MVMR, cis-MR; one-sample vs two-sample bias direction; winner's curse; STROBE-MR |
| colocalization-analysis | Shared-causal-variant inference with coloc.abf / coloc.susie / SMR-HEIDI / eCAVIAR / PWCoCo / moloc / HyPrColoc / SharePro; p12 sensitivity; LD-mismatch diagnostics |
| fine-mapping | Posterior credible-set inference with SuSiE / SuSiE-inf / FINEMAP / PolyFun / SuSiEx; LD-mismatch diagnostics; functional priors; multi-ancestry |
| mediation-analysis | Causal mediation with mediation / CMAverse / HIMA2; 4-way decomposition; MR-mediation / MVMR-mediation; Mediational E-value; high-dimensional EWAS mediators |
| pleiotropy-detection | CHP vs UHP framework; CAUSE / LHC-MR / LCV / MR-Clust / MR-PRESSO / MR-Mix / contamination mixture; Egger NOME / SIMEX; Steiger-filter caveats |
| transcriptome-wide-association | Gene-level TWAS with FUSION / S-PrediXcan / S-MultiXcan / UTMOST / FOCUS / MA-FOCUS; tissue selection; LD-induced false positives; TWAS-MR-coloc triangulation |
| heritability-partitioning | SNP heritability with LDSC / LDAK SumHer / HDL / HESS / BOLT-REML / graphREML; stratified LDSC (Finucane 2015); baseline-LD model; cell-type prioritization |
| proteome-mr-drug-target | Cis-pQTL MR for drug-target validation with UKB-PPP / deCODE / Fenland; Olink vs SomaScan platform discordance; PAV exclusion; phenome-wide on-target adverse effects |
| effector-gene-prioritization | Variant-to-gene mapping with Open Targets L2G / MAGMA / FUMA / cS2G / PoPS / FLAMES; multi-evidence concordance; cross-reference to ABC and ENCODE-rE2G |
| genetic-correlation | Bivariate rg with cross-trait LDSC / HDL / LAVA / rho-HESS / Popcorn; LDSC intercept absorbs sample overlap; local vs global rg; rg as MR validity check |
| genomic-sem | Structural equation modeling of GWAS with GenomicSEM and MTAG; common-factor GWAS with Q_SNP heterogeneity; stratified GenomicSEM; sample overlap via LDSC sampling covariance |

## Example Prompts

- "Test whether BMI causally affects type 2 diabetes using GWAS summary statistics, with a full pleiotropy sensitivity battery"
- "Check if this GWAS signal and an eQTL share a causal variant; if there is allelic heterogeneity, switch to coloc.susie"
- "Run cis-pQTL MR for PCSK9 on coronary artery disease using UKB-PPP, then validate with colocalization and replicate on SomaScan"
- "Fine-map this GWAS locus with susie_rss, include LD-mismatch diagnostics, and apply PolyFun functional priors"
- "Does CpG methylation mediate the effect of SNPs on disease risk across thousands of CpGs? Run HIMA2 with FDR control"
- "Estimate stratified heritability for this trait across baseline-LD annotations and prioritize the most-relevant tissues"
- "Run TWAS-MR-coloc triangulation: S-PrediXcan in liver, cis-eQTL MR, coloc.susie, and FOCUS fine-mapping to nominate the causal gene"
- "Build a common-factor GWAS over depression, anxiety, and PTSD with GenomicSEM and report Q_SNP heterogeneity"
- "Detect correlated horizontal pleiotropy with CAUSE before running standard IVW on this trait pair"
- "Run cross-trait LDSC to estimate genetic correlation between schizophrenia and bipolar disorder, then run LAVA for local rg"

## Requirements

```bash
# R packages from CRAN
install.packages(c('remotes', 'mediation', 'coloc', 'susieR', 'MendelianRandomization',
                   'lavaan', 'HIMA', 'EValue', 'causalweight'))

# GitHub-only R packages
remotes::install_github(c('MRCIEU/TwoSampleMR',
                          'rondolab/MR-PRESSO',
                          'jean997/cause',
                          'LizaDarrous/lhcMR',
                          'cnfoley/mrclust',
                          'qingyuanzhao/mr.raps',
                          'jrs95/hyprcoloc',
                          'jwr-git/pwcoco',
                          'BS1125/CMAverse',
                          'GenomicSEM/GenomicSEM',
                          'zhenin/HDL',
                          'josefin-werme/LAVA',
                          'HDTian/DRMR'))

# Python tooling
pip install metaxcan pyfocus mtag opentargets-genetics
# LDSC python3 fork
pip install git+https://github.com/belowlab/ldsc

# CLI binaries (download from upstream)
# - LDAK 6+ (dougspeed.com)
# - FINEMAP 1.4+ (christianbenner.com)
# - SMR (cnsgenomics.com/software/smr)
# - MAGMA (ctglab.nl/software/magma)
# - SuSiEx (github.com/getian107/SuSiEx)
# - FUSION (gusevlab.org/projects/fusion)
# - BOLT-LMM / BOLT-REML (alkesgroup.broadinstitute.org/BOLT-LMM)
# - GCTA (cnsgenomics.com/software/gcta)
# - HESS (github.com/huwenboshi/hess)
# - Popcorn (github.com/brielin/Popcorn)
```

## Related Skills

- **population-genetics** - GWAS summary statistics, LD, association testing
- **clinical-databases** - ClinVar, gnomAD, dbSNP variant annotation
- **pathway-analysis** - GO enrichment and GSEA for causal gene sets
- **differential-expression** - eQTL data sources (GTEx, eQTLGen) for colocalization, mediation, and TWAS
- **atac-seq** - ABC and ENCODE-rE2G enhancer-gene linking; variant chromatin effect (chromBPNet) for effector-gene prioritization
- **gene-regulatory-networks** - Co-expression context for causal gene nomination
