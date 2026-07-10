---
name: bio-workflows-causal-genomics-pipeline
description: End-to-end post-GWAS causal inference pipeline orchestrating heritability partitioning, genetic correlation, Mendelian randomization with CHP-aware sensitivity (CAUSE / LHC-MR), colocalization, fine-mapping with SuSiE / FOCUS, mediation, TWAS triangulation, cis-pQTL drug-target MR, effector-gene prioritization (L2G / PoPS / cS2G), and GenomicSEM common-factor GWAS. Use when triangulating causal inference across multiple complementary methods, prioritizing tissues via stratified LDSC, nominating or de-risking drug targets, mapping a lead SNP to a candidate effector gene, modeling shared genetic architecture across correlated traits, or producing a STROBE-MR-compliant publication-grade evidence battery from GWAS summary statistics.
tool_type: r
primary_tool: TwoSampleMR
workflow: true
depends_on:
  - causal-genomics/mendelian-randomization
  - causal-genomics/colocalization-analysis
  - causal-genomics/fine-mapping
  - causal-genomics/pleiotropy-detection
  - causal-genomics/mediation-analysis
  - causal-genomics/transcriptome-wide-association
  - causal-genomics/heritability-partitioning
  - causal-genomics/proteome-mr-drug-target
  - causal-genomics/effector-gene-prioritization
  - causal-genomics/genetic-correlation
  - causal-genomics/genomic-sem
qc_checkpoints:
  - after_h2_ldsc: "Mean chi-squared > 1.02; h2 SE < 0.02; intercept ratio < 0.3"
  - after_rg_check: "If abs(rg) > 0.3 then CHP suspected and CAUSE/LHC-MR required"
  - after_instrument_selection: "F-statistic > 10 (two-sample) or > 20 (one-sample); no palindromic SNPs at MAF near 0.5"
  - after_mr: "IVW + Egger + weighted median + weighted mode concordance"
  - after_sensitivity: "MR-PRESSO global p, Egger intercept p (with Isq >= 0.9 for NOME), Steiger directionality"
  - after_chp_check: "CAUSE delta_elpd z > 1.96 OR LHC-MR posterior excludes null"
  - after_coloc: "PP.H4 >= 0.7 triangulation, >= 0.8 publication, >= 0.95 industry; p12 sensitivity stable"
  - after_finemapping: "Credible-set purity (min_abs_corr) >= 0.5; estimate_s_rss lambda < 0.05 if external LD"
  - after_twas: "FOCUS PIP >= 0.8 for candidate causal gene; tissue selected via stratified LDSC"
  - after_effector_gene: "L2G + PoPS + coloc + TWAS concordance >= 3 of 6 evidence streams"
  - after_mediation: "rho_crit > 0.3 OR mediational E-value > 2 (Imai sensitivity)"
---

## Version Compatibility

Reference examples tested with: TwoSampleMR 0.5+, MR-PRESSO 1.0+, coloc 5.2+, susieR 0.12+, MendelianRandomization 0.9+, ldsc 1.0.1 (python3 fork), MetaXcan 0.7+, pyfocus 0.6+, MAGMA 1.10+, MRlap 0.0.3+, cause 1.2+, lhcMR 0.0.0.9000+, HDL 1.4+, LAVA 0.1+, GenomicSEM 0.0.5+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python: `pip show <pkg>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Causal Genomics Pipeline

**"Run post-GWAS causal inference from summary statistics"** -> Orchestrate heritability partitioning and tissue prioritization, genetic-correlation diagnostics, instrument selection, Mendelian randomization with CHP-aware sensitivity, colocalization, fine-mapping (SuSiE / FOCUS), mediation, TWAS triangulation, cis-pQTL drug-target MR, effector-gene prioritization, and (optionally) GenomicSEM common-factor GWAS to triangulate causal evidence and nominate publication-grade causal exposures, genes, and mechanisms.

This is a workflow skill: it owns the chaining decisions and hand-offs, not the internals of any one step. Every step below cross-references the component skill that teaches its mechanism.

## The governing principle

A causal claim is decided at the seams where summary statistics meet, not inside any single method.

1. **Everything shares ONE genome build, and effect alleles are harmonized once.** The exposure and outcome sumstats, the LD reference, and any eQTL/pQTL panel must share a build; if they differ, liftover ONCE, strand-aware (the BBIS inverted-region danger), BEFORE harmonization. `harmonise_data` aligns effect alleles; a flipped palindrome (MAF near 0.5) flips the causal-effect SIGN silently — drop intermediate-frequency palindromes.
2. **The LD reference ancestry must match the GWAS ancestry, committed once and inherited by clumping, coloc.susie, fine-mapping (SuSiE-rss), and TWAS.** A mismatched LD panel corrupts credible sets and colocalization silently; the in-pipeline alarms are `estimate_s_rss` lambda <0.05 and credible-set purity `min_abs_corr` >=0.5.
3. **Order is causal: a diagnostic at step 0 decides whether a later step runs.** Cross-trait LDSC `abs(rg)>0.3` makes correlated horizontal pleiotropy (CHP) suspected, so CAUSE/LHC-MR becomes MANDATORY — MR-PRESSO is BLIND to CHP. Tissue picked at step 0 (stratified LDSC) flows into TWAS; picking it post-hoc is circular. Steiger pre-filter directionality and LD-clump BEFORE coloc/fine-mapping.
4. **Triangulate; refuse the single number.** No single MR estimator, coloc PP.H4, or evidence stream is decisive — report IVW+Egger+median+mode concordance, coloc PP.H4 with a p12 sweep, and effector-gene evidence across >=3 of 6 streams. A bare headline statistic hides the seam.

## Pipeline Overview

```
GWAS Summary Statistics (exposure + outcome)
    |
    v
[0. Pre-flight: h2 + tissue prioritization + rg diagnostic]
    LDSC / S-LDSC baseline-LD / Finucane 2018 cell-type
    Cross-trait LDSC / HDL / LAVA --> if abs(rg) > 0.3 then CHP-aware MR required
    |
    v
[1. Instrument Selection] -----> LD clumping, F-stat filtering, Steiger pre-filter
    |
    v
[2. Mendelian Randomization] --> IVW, MR-Egger, Weighted Median/Mode, MR-RAPS
    |
    +--> [3. Sensitivity] -------> MR-PRESSO, Egger intercept (Isq), leave-one-out, Steiger
    |
    +--> [3b. CHP-aware MR] -----> CAUSE (delta_elpd), LHC-MR posterior (if rg > 0.3)
    |
    v
[4. Colocalization] -----------> coloc.abf / coloc.susie / HyPrColoc / SMR-HEIDI
    |
    v
[5. Fine-Mapping] -------------> SuSiE rss + estimate_s_rss / FINEMAP-inf / PolyFun
    |
    v
[6. Mediation Analysis] -------> Network MR / MVMR / CMAverse 4-way
    |
    +--> [7. TWAS triangulation] -> FUSION / S-PrediXcan / FOCUS PIP >= 0.8
    |
    +--> [8. Cis-pQTL drug-target MR] -> UKB-PPP / deCODE + cross-platform replication
    |
    v
[9. Effector-gene prioritization] -> Open Targets L2G + PoPS + cS2G + coloc + TWAS (>= 3 of 6 evidence)
    |
    v
[10. (optional) GenomicSEM common-factor GWAS] -> factor model + Q_SNP
    |
    v
Triangulated causal-evidence summary across methods
```

## Step 0: Pre-flight - Heritability, Tissue Prioritization, Genetic Correlation

```bash
ldsc.py --h2 trait.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out trait.h2
ldsc.py --h2 trait.sumstats.gz --ref-ld-chr baselineLD.,cts_v8. --w-ld-chr weights. --ref-ld-chr-cts Multi_tissue_gene_expr.ldcts --out trait.cts
ldsc.py --rg trait1.sumstats.gz,trait2.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out rg
```

**Goal:** Confirm heritable signal, pick the right tissue for TWAS / V2G, and detect shared heritable confounding that mandates CHP-aware MR (CAUSE / LHC-MR).

**Reconciliation:** S-LDSC mean chi-squared > 1.02 with h2 SE < 0.02 and intercept ratio < 0.3 is required. Cell-type prioritization with coefficient_p < 0.05 / N_tissues nominates the tissue for downstream TWAS weights and ABC enhancer-gene priors. If cross-trait LDSC abs(rg) > 0.3 (and HDL sample-overlap < 5%), Step 3b becomes mandatory. See causal-genomics/heritability-partitioning and causal-genomics/genetic-correlation.

## Step 1: Instrument Selection

```r
library(TwoSampleMR)
exposure_dat <- read_exposure_data(filename = 'exposure_gwas.tsv', sep = '\t',
    snp_col = 'SNP', beta_col = 'BETA', se_col = 'SE',
    effect_allele_col = 'A1', other_allele_col = 'A2',
    eaf_col = 'EAF', pval_col = 'P')
exposure_dat <- subset(exposure_dat, pval.exposure < 5e-8)
exposure_dat <- clump_data(exposure_dat, clump_r2 = 0.001, clump_kb = 10000)
exposure_dat$F_stat <- (exposure_dat$beta.exposure / exposure_dat$se.exposure)^2
exposure_dat <- subset(exposure_dat, F_stat >= 10)
exposure_dat <- subset(exposure_dat, !(eaf.exposure > 0.42 & eaf.exposure < 0.58 & substr(effect_allele.exposure,1,1) %in% c('A','T') & substr(other_allele.exposure,1,1) %in% c('A','T')))
```

For cis-MR (drug target) use clump_r2 = 0.1 within +/- 500 kb of the gene. Use 5e-9 if M > 5M variants tested.

## Step 2: Mendelian Randomization

```r
outcome_dat <- read_outcome_data(filename = 'outcome_gwas.tsv', sep = '\t',
    snp_col = 'SNP', beta_col = 'BETA', se_col = 'SE',
    effect_allele_col = 'A1', other_allele_col = 'A2',
    eaf_col = 'EAF', pval_col = 'P')
dat <- harmonise_data(exposure_dat, outcome_dat)
mr_results <- mr(dat, method_list = c('mr_ivw', 'mr_egger_regression',
    'mr_weighted_median', 'mr_weighted_mode'))
```

Concordance across IVW, Egger, weighted median, and weighted mode is the headline causal claim. See causal-genomics/mendelian-randomization.

## Step 3: Sensitivity Analysis

```r
library(MRPRESSO)
presso <- mr_presso(BetaOutcome = 'beta.outcome', BetaExposure = 'beta.exposure',
    SdOutcome = 'se.outcome', SdExposure = 'se.exposure',
    OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat,
    NbDistribution = 5000, SignifThreshold = 0.05)
egger_int <- mr_pleiotropy_test(dat)
isq <- Isq(abs(dat$beta.exposure), dat$se.exposure)   # I^2_GX needs same-sign effects; pass abs(beta)
het <- mr_heterogeneity(dat)
loo <- mr_leaveoneout(dat)
steiger <- directionality_test(dat)
```

Isq >= 0.9 is required for the MR-Egger NOME assumption; below that, run SIMEX correction or drop Egger. See causal-genomics/pleiotropy-detection.

## Step 3b: CHP-aware MR (when rg > 0.3 or shared confounder suspected)

```r
library(cause)
library(lhcMR)
cause_fit <- cause(X = cause_dat, variants = top_vars, param_ests = params)
elpd <- summary(cause_fit)$elpd
# lhcMR is a 3-call chain: merge_sumstats (takes LD/rho paths) -> calculate_SP -> lhc_mr
lhc_df <- merge_sumstats(input.files, trait.names, LD.filepath = ld_path, rho.filepath = rho_path)
SP_list <- calculate_SP(lhc_df, trait.names, nStep = 2, SP_single = 3, SP_pair = 50)
lhc_fit <- lhc_mr(SP_list, trait.names, paral_method = 'lapply', nBlock = 200)
```

CAUSE reports delta_elpd of sharing-vs-causal model; z > 1.96 favors true causation over CHP. LHC-MR jointly estimates causal effect and confounder effect via likelihood; the 95% credible interval excluding zero is the causal-effect verdict. Required when LDSC abs(rg) > 0.3. See causal-genomics/pleiotropy-detection.

## Step 4: Colocalization

```r
library(coloc)
d1 <- list(beta = exposure_locus$BETA, varbeta = exposure_locus$SE^2,
    snp = exposure_locus$SNP, position = exposure_locus$BP,
    type = 'quant', N = exposure_n, MAF = exposure_locus$EAF)
d2 <- list(beta = outcome_locus$BETA, varbeta = outcome_locus$SE^2,
    snp = outcome_locus$SNP, position = outcome_locus$BP,
    type = 'cc', N = outcome_n, s = case_fraction, MAF = outcome_locus$EAF)
result <- coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
sens <- coloc::sensitivity(result, rule = 'H4 > 0.7')
```

PP.H4 >= 0.7 for triangulation, >= 0.8 for publication, >= 0.95 for industry-grade target packages. Always sweep p12 over 1e-6 to 5e-5; conclusions must be stable across the sweep. For allelic heterogeneity use coloc.susie. See causal-genomics/colocalization-analysis.

## Step 5: Fine-Mapping with SuSiE

```r
library(susieR)
R <- as.matrix(read.csv('ld_matrix.csv', row.names = 1))
diag_s <- estimate_s_rss(z = locus_stats$BETA / locus_stats$SE, R = R, n = sample_size)
fitted <- susie_rss(bhat = locus_stats$BETA, shat = locus_stats$SE,
    R = R, n = sample_size, L = 10, coverage = 0.95, min_abs_corr = 0.5)
cs <- fitted$sets$cs
```

If `estimate_s_rss` lambda > 0.05 the external LD is mismatched; rerun with in-sample LD or use SuSiE-inf / FINEMAP-inf. min_abs_corr >= 0.5 (r-squared >= 0.25) is the purity threshold for retaining a credible set. For HLA use L = 20-30. See causal-genomics/fine-mapping.

## Step 6: Mediation Analysis

```r
library(TwoSampleMR)
mv_exposures <- mv_extract_exposures(c('ieu-a-2', 'ieu-a-1089'))
mv_outcome <- extract_outcome_data(mv_exposures$SNP, 'ieu-a-7')
mvdat <- mv_harmonise_data(mv_exposures, mv_outcome)
mvmr_result <- mv_multiple(mvdat)
```

Indirect effect = total - direct. For molecular mediators (expression, methylation, protein), prefer two-step MR with cis-instruments at the mediator. Run Imai sensitivity (rho_crit) or mediational E-value > 2. See causal-genomics/mediation-analysis.

## Step 7: TWAS Triangulation

```bash
python MetaXcan/SPrediXcan.py --model_db_path mashr_Whole_Blood.db \
    --covariance mashr_Whole_Blood.txt.gz --gwas_file gwas.txt.gz \
    --snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 \
    --beta_column BETA --se_column SE --pvalue_column P \
    --output_file twas.csv
focus finemap gwas.sumstats.gz 1000G.EUR.QC.1 mashr.db --chr 1 --locations 37:EUR --out twas.focus
```

**Goal:** Nominate gene-level causal hits and prune LD-induced TWAS false positives.

Tissue is picked from Step 0 stratified LDSC; Bonferroni at 0.05 / N_tissues. FOCUS PIP >= 0.8 retains a single candidate causal gene per region; without FOCUS, co-regulated TWAS hits cannot be distinguished. Cross-reference TWAS hits with coloc.susie PP.H4 and cis-eQTL MR for triangulation. See causal-genomics/transcriptome-wide-association.

## Step 8: Cis-pQTL Drug-Target MR

```r
library(TwoSampleMR)
pqtl_dat <- extract_instruments(outcomes = 'prot-a-XXX', p1 = 5e-8, clump = TRUE,
    r2 = 0.1, kb = 1000)
pqtl_dat <- subset(pqtl_dat, chr.exposure == target_chr & pos.exposure > target_tss - 500000 & pos.exposure < target_tss + 500000)   # extract_instruments returns chr.exposure/pos.exposure
out_dat <- extract_outcome_data(snps = pqtl_dat$SNP, outcomes = c('ieu-a-7', 'ieu-b-31'))
dat <- harmonise_data(pqtl_dat, out_dat)
mr_results <- mr(dat, method_list = c('mr_ivw', 'mr_wald_ratio'))
```

**Goal:** Mimic pharmacological inhibition of a drug target via cis-pQTL and triangulate with coloc.

Cross-platform replication on Olink (UKB-PPP) and SomaScan (deCODE) is mandatory; the two platforms are concordant for ~60% of proteins and discordant calls are platform artifacts. Run pheWAS for on-target adverse effects, PAV-excluded sensitivity, and coloc.susie PP.H4 >= 0.8 at the cis-pQTL locus. See causal-genomics/proteome-mr-drug-target.

## Step 9: Effector-Gene Prioritization

```bash
magma --bfile g1000_eur --pval gwas.tsv N=N --gene-annot genes.annot --out trait
python munge_feature_directory.py --gene_annot_path genes.txt --feature_dir features/ --save_prefix pops
python pops.py --gene_annot_path genes.txt --feature_mat_prefix pops --num_feature_chunks 2 --magma_prefix trait --out_prefix trait.pops
```

**Goal:** Map each fine-mapped credible set to a candidate effector gene by integrating six evidence streams.

Integrate: (1) Open Targets L2G (Mountjoy 2021), (2) PoPS similarity score (Weeks 2023), (3) cS2G combined SNP-to-gene (Gazal 2022), (4) coloc.susie PP.H4 with eQTL/pQTL, (5) FOCUS TWAS PIP, (6) ABC / ENCODE-rE2G enhancer-gene linking. Require >= 3 of 6 concordant evidence streams for high-confidence claim. L2G and PoPS disagree by design (different feature regimes); report both. See causal-genomics/effector-gene-prioritization.

## Step 10 (optional): GenomicSEM Common-Factor GWAS

```r
library(GenomicSEM)
ldsc_output <- ldsc(traits = c('t1.sumstats.gz','t2.sumstats.gz','t3.sumstats.gz'),
    sample.prev = c(NA, NA, NA), population.prev = c(NA, NA, NA),
    ld = ld_path, wld = ld_path, trait.names = c('t1','t2','t3'))
model <- 'F1 =~ NA*t1 + t2 + t3\nF1 ~~ 1*F1'
fit <- usermodel(ldsc_output, model = model, estimation = 'DWLS')
factor_gwas <- userGWAS(covstruc = ldsc_output, SNPs = sumstats_combined,
    model = paste0(model, '\nF1 ~ SNP\nt1 + t2 + t3 ~ 0*SNP'),
    estimation = 'DWLS', sub = c('F1~SNP'))
```

Heywood cases (negative residual variance) require fixing residuals positive or dropping the indicator. Verify CFI > 0.95 and RMSEA < 0.06. Q_SNP p > 0.05 confirms factor-level (not trait-specific) signal. See causal-genomics/genomic-sem.

## Parameter Recommendations

| Step | Parameter | Recommendation |
|------|-----------|----------------|
| Instruments | p-value | 5e-8 standard; 5e-9 if M > 5M variants tested |
| Instruments | F-statistic | >= 10 two-sample; >= 20 one-sample |
| Instruments | clump_r2 | 0.001 polygenic; 0.1 cis-MR |
| Instruments | clump_kb | 10000 (10 Mb) |
| Coloc | p12 prior | 1e-5 standard; 5e-6 conservative; 1e-6 trans-eQTL |
| Coloc | PP.H4 | >= 0.7 triangulation; >= 0.8 publication; >= 0.95 industry |
| SuSiE | L | 10 default; 20-30 HLA |
| SuSiE | coverage | 0.95 standard; 0.9 if N < 1000 |
| SuSiE | min_abs_corr | 0.5 default (r-squared >= 0.25) |
| MR-PRESSO | NbDistribution | 1000 exploratory; >= 5000 publication; >= 10000 stringent |
| TWAS | FOCUS PIP | >= 0.8 candidate causal gene |
| TWAS | tissue Bonferroni | 0.05 / N_tissues (~2.5e-4 for 200 tissues) |
| LDSC | mean chi-squared | > 1.02 for h2 interpretability |
| LDSC | rg trigger for CHP-MR | abs(rg) > 0.3 |
| LDSC | HDL sample overlap | < 5% |

## Common Errors

| Issue | Likely Cause | Solution |
|-------|--------------|----------|
| Nonsense / sign-flipped MR or coloc | Build/liftover mismatch across exposure/outcome/LD/eQTL | One build; liftover once strand-aware (BBIS danger), THEN harmonize |
| Causal-effect sign flipped | Strand-flip/allele-swap at a MAF~0.5 palindrome in harmonise | Drop intermediate-frequency palindromes; verify effect-allele alignment |
| Corrupted credible sets / wrong coloc | LD-panel ancestry mismatch (clumping/coloc/fine-map/TWAS) | Ancestry-matched LD; check estimate_s_rss lambda; prefer in-sample LD |
| "Causal" effect passes sensitivity but is confounded | MR-PRESSO blind to CHP; rg>0.3 not checked | Run the rg gate at step 0; CAUSE/LHC-MR when triggered |
| Two-sample MR biased toward observational | Sample overlap between exposure and outcome GWAS | Check overlap (<5% via HDL); MRlap/HDL correction |
| No instruments | Underpowered GWAS | Relax p-value to 5e-6 with caution |
| Weak instruments (F < 10) | Small effect SNPs | Drop weak instruments; use better-powered GWAS |
| Inconsistent MR methods | Pleiotropy | Check MR-PRESSO outliers; use weighted median |
| Egger intercept p < 0.05 | Directional pleiotropy | Report Egger estimate; check Isq >= 0.9 |
| PP.H3 > PP.H4 | Different causal variants | Use coloc.susie for allelic heterogeneity |
| No credible sets | LD matrix issues | Check estimate_s_rss lambda; use in-sample LD |
| Steiger reverse direction | Reverse causation | Run bidirectional MR; pre-filter on Steiger |
| MR-PRESSO blind to apparent confounder | Correlated horizontal pleiotropy | Run CAUSE or LHC-MR (see pleiotropy-detection) |
| TWAS hit at gene-dense locus | LD-induced false positive | Run FOCUS fine-mapping; require PIP >= 0.8 |
| Cis-pQTL MR positive, replication fails | Olink/SomaScan platform discordance | Cross-platform replication; PAV-excluded sensitivity |
| L2G + PoPS disagree | Different feature regimes | Report both; require concordance for high-confidence claim |
| Q_SNP heterogeneity in common-factor GWAS | Factor mis-specification | userGWAS with per-trait paths; report Q_pval |

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
- workflows/gwas-pipeline - Upstream: produces the sumstats (CHR/POS/EA/OA/EAF/BETA/SE/P/N, build documented) this pipeline consumes

## References

- Hemani G, Zheng J, Elsworth B, et al (2018) The MR-Base platform supports systematic causal inference across the human phenome. *eLife* 7:e34408. DOI 10.7554/eLife.34408. (TwoSampleMR / harmonisation.)
- Giambartolomei C, Vukcevic D, Schadt EE, et al (2014) Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. *PLoS Genetics* 10:e1004383. DOI 10.1371/journal.pgen.1004383. (coloc.)
- Wang G, Sarkar A, Carbonetto P, Stephens M (2020) A simple new approach to variable selection in regression, with application to genetic fine mapping. *Journal of the Royal Statistical Society Series B* 82:1273-1300. DOI 10.1111/rssb.12388. (SuSiE.)
- Bulik-Sullivan BK, Loh PR, Finucane HK, et al (2015) LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. *Nature Genetics* 47:291-295. DOI 10.1038/ng.3211. (LDSC intercept.)
