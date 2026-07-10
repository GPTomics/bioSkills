# Causal Genomics Pipeline - Usage Guide

## Overview
Complete post-GWAS causal inference workflow that triangulates evidence across heritability partitioning, genetic correlation, Mendelian randomization (with CHP-aware sensitivity via CAUSE / LHC-MR), colocalization, fine-mapping (SuSiE / FOCUS), mediation, TWAS, cis-pQTL drug-target MR, effector-gene prioritization (L2G / PoPS / cS2G), and (optionally) GenomicSEM common-factor GWAS. Suitable for nominating causal exposures, picking the right tissue for downstream analyses, mapping a lead SNP to a candidate effector gene, de-risking or nominating a drug target, and producing a STROBE-MR-compliant publication-grade evidence battery from GWAS summary statistics.

## Prerequisites
```r
install.packages('remotes')
remotes::install_github('MRCIEU/TwoSampleMR')
remotes::install_github('rondolab/MR-PRESSO')
remotes::install_github('jean997/cause')
remotes::install_github('zhenin/HDL/HDL')
remotes::install_github('josefin-werme/LAVA')
remotes::install_github('jrs95/MRlap')
remotes::install_github('LizaDarrous/lhcMR')
install.packages('coloc')
install.packages('susieR')
install.packages('MendelianRandomization')
remotes::install_github('GenomicSEM/GenomicSEM')   # GenomicSEM is GitHub-only (not on CRAN)
```

```bash
git clone https://github.com/bulik/ldsc.git
conda env create -n ldsc -f ldsc/environment.yml
git clone https://github.com/hakyimlab/MetaXcan.git
pip install pyfocus
# MAGMA ships as a zip, not a git repo. Grab the current program zip from the official CNCR page
# (https://cncr.nl/research/magma/) -- the old vu.data.surfsara.nl direct links now 404 after SURF's
# domain migration, so fetch the link from that page rather than hard-coding one here.
unzip magma_v1.10.zip
git clone https://github.com/FinucaneLab/pops.git   # PoPS is run as pops.py from the clone; PyPI `pops` is unrelated
```

**Input data:**
- GWAS summary statistics for exposure and outcome (TSV with SNP, BETA, SE, A1, A2, EAF, P, N columns)
- Pre-computed LDSC weights and baseline-LD annotation (eur_w_ld_chr, baselineLD)
- TWAS prediction weights (PredictDB MASHR models or FUSION elastic net) matched to the tissue prioritized in Step 0
- LD matrix for fine-mapping loci (in-sample preferred; reference-panel acceptable with estimate_s_rss check)
- Cis-pQTL summary stats for drug-target MR (UKB-PPP, deCODE, Fenland, INTERVAL)
- Open Targets L2G features and PoPS gene-feature table for effector-gene step

## Quick Start
Tell your AI agent what you want to do:
- "Run the full post-GWAS causal-inference pipeline on my GWAS summary stats"
- "Estimate heritability, prioritize tissues, then run MR with CHP-aware sensitivity"
- "Triangulate TWAS, colocalization, and cis-eQTL MR to nominate the causal gene"
- "Run cis-pQTL drug-target MR with Olink and SomaScan cross-platform replication"
- "Prioritize an effector gene from a GWAS lead locus using L2G, PoPS, and ABC"
- "Test a common-factor GWAS model across my correlated traits via GenomicSEM"

## Example Prompts

### Full Pipeline
> "I have GWAS summary statistics for LDL cholesterol and coronary heart disease. Run the complete causal-inference pipeline including heritability partitioning, tissue prioritization, MR with CHP-aware sensitivity, coloc, fine-mapping, and effector-gene prioritization."

> "Perform MR, sensitivity analysis, and colocalization to test if BMI causally affects type 2 diabetes, with CAUSE if rg > 0.3."

### Heritability and Tissue Prioritization
> "Estimate SNP heritability with LDSC, partition with the baseline-LD model, and prioritize tissues via Finucane 2018 cell-type LDSC for my schizophrenia GWAS."

> "Run stratified LDSC on my GWAS to pick the right GTEx tissue for downstream TWAS."

### Genetic Correlation
> "Compute cross-trait LDSC rg between my exposure and outcome to decide whether CAUSE or LHC-MR is required."

> "Run LAVA local genetic correlation across the genome to find regions of shared signal."

### MR and Sensitivity
> "Select strong instruments and check F-statistics for my exposure GWAS."

> "Run MR-PRESSO with NbDistribution = 10000 to detect outlier instruments."

> "Run CAUSE because LDSC rg > 0.3 between my traits and CHP is plausible."

### Colocalization and Fine-Mapping
> "Colocalize my GWAS locus with whole-blood eQTL using coloc.susie because the region has multiple signals."

> "Fine-map my GWAS locus with SuSiE rss, run estimate_s_rss for LD diagnostics, and report 95% credible sets with min_abs_corr >= 0.5."

### TWAS Triangulation
> "Run S-PrediXcan in whole-blood MASHR models and fine-map TWAS hits with FOCUS at PIP >= 0.8."

> "Triangulate TWAS, coloc.susie, and cis-eQTL MR to nominate a single causal gene per locus."

### Drug-Target MR
> "Run cis-pQTL MR for IL6R using UKB-PPP and replicate in deCODE SomaScan with PAV-excluded sensitivity."

> "Screen on-target adverse phenotypes with pheWAS at the cis-pQTL for my drug target."

### Effector-Gene Prioritization
> "Map this GWAS lead variant to a candidate effector gene using Open Targets L2G, PoPS, cS2G, coloc, and ABC."

> "Reconcile L2G and PoPS calls at my locus and report concordance across 6 evidence streams."

### Mediation
> "Test if CRP mediates the causal effect of BMI on CHD using two-step MR or MVMR with Imai sensitivity."

> "Run CMAverse 4-way decomposition for exposure-mediator interaction."

### GenomicSEM
> "Fit a common-factor GenomicSEM model across my 5 psychiatric GWAS and run a common-factor GWAS with Q_SNP heterogeneity testing."

> "Compare GenomicSEM common-factor GWAS against MTAG for the same trait set."

## What the Agent Will Do
1. Estimate SNP heritability with LDSC and run partitioned S-LDSC for functional categories
2. Prioritize tissues via Finucane 2018 cell-type-specific LDSC (informs Step 7 TWAS tissue choice)
3. Compute cross-trait LDSC rg; if abs(rg) > 0.3, mandate CHP-aware MR (CAUSE / LHC-MR) in Step 3b
4. Select genome-wide-significant independent instruments with LD clumping, F-statistic, palindromic filters, and Steiger pre-filter
5. Run primary MR across IVW, MR-Egger, weighted median, and weighted mode for concordance
6. Run sensitivity battery: MR-PRESSO, Egger intercept (with Isq >= 0.9 for NOME), leave-one-out, Cochran Q, Steiger directionality
7. If rg > 0.3, run CAUSE delta_elpd and LHC-MR posterior to distinguish causation from CHP
8. Colocalize the locus with coloc.abf or coloc.susie; sweep p12 prior; require PP.H4 >= 0.7 (triangulation) or 0.8 (publication)
9. Fine-map with susie_rss; run estimate_s_rss for LD diagnostics; require min_abs_corr >= 0.5
10. Mediation: two-step MR or MVMR with Imai rho_crit / mediational E-value sensitivity
11. TWAS: run S-PrediXcan in the prioritized tissue, fine-map with FOCUS PIP >= 0.8
12. Cis-pQTL drug-target MR: replicate across Olink and SomaScan platforms; run pheWAS and PAV sensitivity
13. Effector-gene prioritization: integrate L2G, PoPS, cS2G, coloc, TWAS, and ABC; require >= 3 of 6 concordant streams
14. Optional GenomicSEM: fit common-factor model, verify CFI > 0.95 and RMSEA < 0.06, run userGWAS with Q_SNP
15. Triangulate evidence and produce STROBE-MR-compliant summary

## Tips
- Always use multiple MR methods; if directions disagree, investigate pleiotropy and run CAUSE / LHC-MR
- F-statistic >= 10 is minimum for two-sample MR; use >= 20 for one-sample MR to avoid weak-instrument bias toward the observational estimate
- Filter palindromic SNPs with MAF near 0.5 (0.42 to 0.58); they cannot be reliably harmonised
- MR-PRESSO global test p > 0.05 means no significant horizontal pleiotropy; below that, report the outlier-corrected estimate, not the raw IVW
- Isq < 0.9 invalidates the MR-Egger NOME assumption; run SIMEX or drop Egger from the headline
- Cross-trait LDSC abs(rg) > 0.3 triggers CHP-aware sensitivity (CAUSE delta_elpd z > 1.96 or LHC-MR posterior excluding zero); HDL is biased when sample overlap is above 5%
- PP.H4 >= 0.7 is the triangulation threshold, >= 0.8 publication-grade, >= 0.95 industry-grade
- Always sweep coloc.abf p12 over 1e-6 to 5e-5; conclusions must be stable across the sweep
- Use coloc.susie (not coloc.abf) for loci with allelic heterogeneity (multiple independent signals)
- Fine-mapping requires accurate LD; check estimate_s_rss lambda < 0.05 before trusting external reference-panel LD
- min_abs_corr >= 0.5 (r-squared >= 0.25) is the SuSiE purity threshold; lower sets are unreliable
- For TWAS, pick the tissue from Step 0 stratified LDSC; do not run all 49 GTEx tissues by default
- FOCUS PIP >= 0.8 is required to retain a single candidate causal gene per region; without it, co-regulated TWAS hits cannot be distinguished
- Cis-pQTL MR requires cross-platform replication on Olink (UKB-PPP) and SomaScan (deCODE); ~40% of proteins are platform-discordant
- For effector-gene prioritization, require >= 3 of 6 concordant evidence streams (L2G, PoPS, cS2G, coloc, TWAS, ABC); L2G and PoPS disagree by design and both should be reported
- Steiger pre-filtering of instruments (variance explained in exposure > variance explained in outcome) reduces reverse-causation bias
- GenomicSEM common-factor GWAS requires CFI > 0.95 and RMSEA < 0.06; Q_SNP p > 0.05 confirms factor-level (not trait-specific) signal
- For mediation, run Imai sensitivity (rho_crit > 0.3 robust) or mediational E-value > 2; molecular mediators are best tested via two-step MR with cis-instruments

## Related Skills
- causal-genomics/mendelian-randomization - IVW, Egger, MR-RAPS, MVMR
- causal-genomics/colocalization-analysis - coloc.abf, coloc.susie, HyPrColoc, SMR-HEIDI
- causal-genomics/fine-mapping - SuSiE rss, FINEMAP-inf, PolyFun, SuSiEx
- causal-genomics/pleiotropy-detection - MR-PRESSO, CAUSE, LHC-MR, contamination-mixture
- causal-genomics/mediation-analysis - Two-step MR, MVMR, CMAverse 4-way, HIMA
- causal-genomics/transcriptome-wide-association - FUSION, S-PrediXcan, FOCUS, UTMOST
- causal-genomics/heritability-partitioning - LDSC, S-LDSC, LDAK, HDL, HESS
- causal-genomics/proteome-mr-drug-target - UKB-PPP, deCODE, cis-pQTL MR, pheWAS
- causal-genomics/effector-gene-prioritization - L2G, PoPS, cS2G, MAGMA, FLAMES
- causal-genomics/genetic-correlation - Cross-trait LDSC, HDL, LAVA, Popcorn
- causal-genomics/genomic-sem - Common-factor GWAS, Q_SNP, MTAG reconciliation
- population-genetics/association-testing - Upstream GWAS methods
- atac-seq/enhancer-gene-linking - ABC / ENCODE-rE2G priors for effector-gene step
- single-cell/preprocessing - scRNA / scATAC tissue priors for stratified LDSC
