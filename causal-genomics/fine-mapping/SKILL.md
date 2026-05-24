---
name: bio-causal-genomics-fine-mapping
description: Resolves GWAS associations to candidate causal variants and credible sets via SuSiE, susie_rss, FINEMAP, CAVIAR, DAP-G, PAINTOR, PolyFun, SuSiEx, MultiSuSiE, and FOCUS. Use when narrowing a GWAS lead SNP to a 95 percent credible set, choosing between in-sample and reference LD, calibrating non-sparse loci with SuSiE-inf or FINEMAP-inf, integrating functional priors via PolyFun, fine-mapping across ancestries with SuSiEx, diagnosing LD mismatch via estimate_s_rss and kriging_rss, handling HLA or long-range LD, or feeding credible sets into coloc.susie for colocalization.
tool_type: r
primary_tool: susieR
---

## Version Compatibility

Reference examples tested with: susieR 0.12.27+, coloc 5.2.3+, FINEMAP 1.4.2+, PolyFun (head of `omerwe/polyfun` 2024), PAINTOR V3.0, SuSiEx (head of `getian107/SuSiEx`), DAP-G (head of `xqwen/dap`), pyfocus 0.8+, R 4.3+, PLINK 1.9 / 2.0.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('susieR')` then `?susie_rss` to confirm argument names (e.g., `prior_weights` vs `prior_variance` semantics)
- CLI: `finemap --help`, `SuSiEx --help`, `PAINTOR --help`, `dap-g --help` to confirm flags
- Python: `polyfun.py --help`

If a call throws an error about an argument that no longer exists, introspect the installed function and adapt rather than retrying.

# Fine-Mapping

**"Narrow my GWAS locus to the variants likely to be causal"** -> Fit a sparse Bayesian regression that propagates LD into posterior inclusion probabilities (PIPs) and credible sets, then validate that credible sets correspond to physically reasonable haplotypes given the LD reference.

- R (summary statistics + LD): `susieR::susie_rss(z, R, n, L=10)` + `estimate_s_rss` LD diagnostic
- R (individual-level genotypes): `susieR::susie(X, y, L=10)`
- CLI (shotgun stochastic search): `finemap --sss --in-files master.z --n-causal-snps 5 --sss-tol 0.001`
- CLI (cross-ancestry joint): `SuSiEx --sst_file=eur.sst,eas.sst --n_gwas=N1,N2 --ref_file=eur.bim,eas.bim --ld_file=eur_ld,eas_ld --pop=EUR,EAS --chr_col=1,1 --snp_col=2,2 --bp_col=3,3 --a1_col=4,4 --a2_col=5,5 --eff_col=6,6 --se_col=7,7 --pval_col=8,8 --chr=<chr> --bp=<start,end> --out_dir=<dir> --out_name=<name>` (column-number flags and `--ld_file` are required; see `SuSiEx --help`)
- Python (functional priors): `polyfun.py --compute-h2-L2` -> per-SNP priors -> susie_rss with `prior_weights=`
- Python (TWAS fine-mapping): `focus finemap` on gene-level Z-scores

Fine-mapping is a Bayesian model selection problem; LD is not noise but structured prior information. Most failure modes trace back to one of three issues: (a) LD reference mismatched to the GWAS sample; (b) the sparse-effects prior being wrong for the locus (polygenic background); or (c) too small an L cap. The `estimate_s_rss()` lambda and `kriging_rss()` per-SNP diagnostic catch (a) before downstream credible sets are reported.

## Algorithmic Taxonomy

| Tool | Model | Input | Strength | Fails when |
|------|-------|-------|----------|------------|
| SuSiE / susie_rss (Wang 2020 JRSSB 82:1273; Zou 2022 PLoS Genet) | Iterative Bayesian sum-of-single-effects (IBSS), variational | Individual-level (X, y) or (z, R, n) | Fast; native PIP + credible sets; pluggable priors; default in modern pipelines | Reference LD mismatched to GWAS sample; locus dominated by polygenic background; >L true effects |
| SuSiE-inf / FINEMAP-inf (Cui 2024 Nat Genet 56:162) | SuSiE + infinitesimal random-effect component | (z, R, n) | Calibrated credible sets when locus is non-sparse (polygenic shoulder around a sparse causal); recommended for biobank-scale GWAS | Very small loci with truly sparse architecture (over-conservative); slower convergence |
| FINEMAP (Benner 2016 Bioinformatics 32:1493) | Shotgun stochastic search over causal configurations | .z + .ld + .master files | Exact Bayes factors at small k; widely cited | Slow at L > 5; binary install only (christianbenner.com); same LD-mismatch fragility as SuSiE |
| CAVIAR / CAVIARBF (Hormozdiari 2014 Genetics 198:497) | Exhaustive enumeration up to k causals | (z, R) | Exact posterior at small k | Combinatorial explosion beyond k=6; legacy method largely superseded by SuSiE |
| DAP-G (Wen 2016 AJHG 98:1114) | Deterministic posterior approximation with adaptive scan | SBAMS format; TORUS for enrichment priors | Fast at QTL scale (whole-transcriptome); pairs with TORUS hierarchical priors | SBAMS format is awkward; less ubiquitous tooling |
| PAINTOR (Kichaev 2014 PLoS Genet 10:e1004722) | EM with binary functional annotations | (z, R, A) per locus | Locus-level functional priors; multi-trait variant | Single-trait mode often matched by PolyFun + SuSiE; slower than SuSiE |
| PolyFun + SuSiE/FINEMAP (Weissbrod 2020 Nat Genet 52:1355) | Stratified LDSC genome-wide -> per-SNP prior_weights | GWAS sumstats + pre-baked baseline-LF | Most powerful single-trait functional prior; PIPs sharpen 1.5-3x | Requires matched-ancestry baseline-LF; runs in two stages |
| SuSiEx (Yuan 2024 Nat Genet 56:1841) | Joint cross-ancestry SuSiE; shared causal, population-specific LD | Per-pop sumstats + per-pop LD reference | Smaller credible sets than per-ancestry meta or marginal fine-mapping; principled when causal variants are shared | Trans-ethnic heterogeneity violated (population-specific causals); ancestry must be cleanly assigned |
| MultiSuSiE (Tashman 2024 medRxiv / 2025 Nat Genet) | Cross-ancestry SuSiE variant; flexible heterogeneity | Per-pop sumstats + per-pop LD | Similar to SuSiEx; alternative implementation | Same as SuSiEx; newer, less battle-tested |
| FOCUS / MA-FOCUS (Mancuso 2019 Nat Genet 51:675) | Probabilistic TWAS fine-mapping over gene models | TWAS Z-scores + gene LD (predicted expression) | Identifies likely causal gene among co-regulated TWAS hits; cross-ancestry MA-FOCUS variant | Requires pre-computed expression weights (e.g., FUSION/PrediXcan); gene-level rather than variant-level inference |

Methodology evolves; verify the latest susieR vignette and the SuSiE-inf paper before locking on a single method. Wang Lab maintains susieR; the IBSS algorithm is stable but argument semantics (e.g., `prior_weights` vs `prior_variance`) have changed across versions.

## Decision Tree by Experimental Scenario

| Scenario | Recommended workflow | Why |
|----------|---------------------|-----|
| Individual-level genotypes available (UKB, in-house cohort) | `susie(X, y, L=10)` | In-sample LD is exact; no mismatch fragility |
| Summary statistics only, ancestry matches reference panel | `susie_rss(z, R, n, L=10)` + `estimate_s_rss` diagnostic | Standard external-LD pattern; verify lambda < 0.05 |
| Single-locus EUR GWAS, sparse architecture | susie_rss with L=10, baseline functional priors optional | Most-common setting; SuSiE default works |
| Locus with strong polygenic shoulder (biobank scale) | SuSiE-inf (Cui 2024) | Adds infinitesimal component; calibrates non-sparse PIPs |
| Multi-ancestry GWAS (EUR + EAS + AFR) | SuSiEx with per-pop sumstats and LD | Joint inference shrinks credible sets; per-ancestry meta loses LD information |
| Locus with > 5 expected independent signals (HLA, lipid loci) | susie_rss with L=20-30 | Default L=10 caps signal count; HLA needs extension |
| TWAS hits with co-regulated genes | FOCUS / MA-FOCUS | Variant-level fine-mapping cannot distinguish co-regulated gene candidates |
| Want functional priors (coding, conserved, regulatory) | PolyFun -> susie_rss with `prior_weights` | Genome-wide SLDSC priors sharpen PIPs more than locus-level annotations |
| QTL fine-mapping (eQTL, sQTL, caQTL) at transcriptome scale | DAP-G + TORUS OR susie_rss per gene | DAP-G is built for QTL throughput; SuSiE works per gene |
| Low-N QTL (GTEx tissue panel, N < 1000) | susie_rss with `coverage = 0.9` (or 0.8); document choice | Default 0.95 returns very wide credible sets at low power; report the relaxed coverage explicitly in methods |
| HLA region (chr6:28-34 Mb) or chr8 inversion | Specialized workflow: stratify haplotypes; consider HLA-specific imputation; or exclude | LD structure is too complex; standard methods unreliable |
| Cross-feed into colocalization | susie_rss -> coloc.susie() | Modern coloc operates on credible sets, not single SNPs |

## Critical LD Diagnostic Block (susie_rss)

**Goal:** Detect LD reference mismatch before reporting credible sets.

**Approach:** `estimate_s_rss()` quantifies the global Z-score / LD inconsistency as a scalar; `kriging_rss()` identifies individual SNPs whose Z-scores are inconsistent with the LD reference (typically genotyping errors, strand flips, or wrong reference panel).

```r
library(susieR)
s_hat <- estimate_s_rss(z = z_scores, R = ld_matrix, n = N)
# s_hat is the inferred scale of LD inconsistency.
# Source: susieR vignette "Diagnostic for summary statistic"; Zou 2022 PLoS Genet.
# Rule of thumb: s_hat < 0.05 acceptable; 0.05-0.10 marginal; > 0.10 refit or change LD reference.

cond_z <- kriging_rss(z = z_scores, R = ld_matrix, n = N)
# cond_z$conditional_dist returns per-SNP expected vs observed z; flag |z_obs - z_exp| > 3
# Common cause: strand flip, allele coding mismatch, or single-SNP imputation error.

# If diagnostic fails: refit with explicit scale parameter to absorb LD mismatch
fit <- susie_rss(z = z_scores, R = ld_matrix, n = N, L = 10, estimate_residual_variance = TRUE)
```

Skipping this block is the dominant cause of irreproducible fine-mapping. Always run before reporting credible sets.

## Per-Tool Failure Modes

### LD reference mismatch (most common)

**Trigger:** External LD matrix from 1000 Genomes / UK Biobank reference used for a GWAS conducted on a different cohort or ancestry mix.

**Mechanism:** Z-scores reflect the GWAS sample's LD; the reference R does not. The susie_rss likelihood depends on `z' R^{-1} z` being consistent with the modeled effects, and inconsistency manifests as spurious credible sets containing tag SNPs from the reference but not from the discovery cohort.

**Symptom:** `estimate_s_rss()` lambda > 0.05; `kriging_rss()` flags many SNPs with `|z_obs - z_exp| > 3`; credible sets contain physically distant SNPs (anti-correlated in LD with the lead) or include all SNPs at the locus.

**Fix:** Use in-sample LD whenever the cohort genotypes are accessible (compute with `plink --r2 square` on the GWAS samples themselves). When only summary statistics are available, ancestry-stratify the LD reference exactly (e.g., 1000G EUR FIN+CEU+GBR+IBS+TSI for a Northern European GWAS, not full EUR). For mixed-ancestry GWAS, fine-map per ancestry then meta-analyze, or move to SuSiEx.

### Non-sparse architecture (biobank scale)

**Trigger:** Locus with one strong signal plus hundreds of weakly associated SNPs (polygenic shoulder); typical at biobank scale.

**Mechanism:** Vanilla SuSiE assumes a sparse sum-of-single-effects prior. With polygenic background, the model misallocates effects, producing inflated credible sets or many small spurious ones. Cui 2024 (Nat Genet 56:162) showed PIPs from SuSiE in this regime are systematically miscalibrated.

**Symptom:** Many small credible sets (5-15 per locus); replication in independent cohorts fails for non-lead credible sets; PIP distribution has a heavy tail.

**Fix:** Use SuSiE-inf or FINEMAP-inf (Cui 2024). These augment the sum-of-single-effects with an infinitesimal random-effect component that absorbs polygenic background. Source: github.com/FinucaneLab/fine-mapping-inf.

### L too small

**Trigger:** Locus with > 5 independent signals (HLA region, APOC1/APOE, LPA, IL6R region for some traits).

**Mechanism:** SuSiE assumes at most L independent effects. When the true number exceeds L, some signals are absorbed into existing components, distorting PIPs and credible sets for the captured signals.

**Symptom:** `length(fit$sets$cs)` equals L (all L slots used); credible set purity for higher-indexed sets is low (`fit$sets$purity[,'min.abs.corr'] < 0.5`); fits with larger L change top-PIP variants.

**Fix:** Increase L iteratively (L=10 -> 20 -> 30) until `length(fit$sets$cs)` < L (susieR auto-prunes unsupported effects so the returned CS count is the effective L). For HLA, start at L=30. The cost is mostly computational, not statistical: SuSiE prunes unused slots, so L=30 is safe when L=10 was right.

### prior_weights vs prior_variance confusion (PolyFun integration)

**Trigger:** Passing PolyFun output to susie_rss with `prior_variance=polyfun_priors` (wrong argument).

**Mechanism:** `prior_variance` in susie_rss is a single scalar (or vector of length L) for the per-effect variance, NOT a per-SNP probability. `prior_weights` is the per-SNP causal probability vector (sums to ~1). Passing PolyFun's per-SNP prior to `prior_variance` is silently accepted but applies a numerically nonsensical per-effect variance.

**Symptom:** PIPs nearly identical to the uniform-prior fit; functional annotations appear to have no effect.

**Fix:** Use `prior_weights = polyfun_priors$snpvar`. Verify with `?susie_rss` in the installed version. Reference: github.com/omerwe/polyfun README, Weissbrod 2020 supplementary methods.

### Credible-set misinterpretation

**Trigger:** Reporting per-variant PIP without distinguishing "in credible set" from "high PIP".

**Mechanism:** The 95 percent credible-set guarantee is `P(causal variant in set) >= 0.95`. Per-variant PIPs within a set do not necessarily sum to 1 across all variants, and PIPs across overlapping sets can double-count posterior mass.

**Symptom:** Reporting "the top PIP variant" when the credible set is wide (size > 50); claiming a single variant is causal when the set contains 30 high-LD SNPs.

**Fix:** Always report (a) number of credible sets, (b) size of each set, (c) purity (`fit$sets$purity[,'min.abs.corr']`), (d) the top PIP variant within the set as the candidate lead. The credible set is the unit of inference; the top PIP variant is a candidate, not a conclusion.

### Cross-ancestry with single-ancestry LD

**Trigger:** Multi-ancestry meta-analyzed GWAS, then susie_rss with EUR LD.

**Mechanism:** Meta-analysis z-scores reflect a weighted mix of population LD structures; no single-population LD matrix matches.

**Fix:** Move to SuSiEx (joint cross-ancestry SuSiE; Yuan 2024). Per-ancestry fine-mapping followed by manual merging loses the shared-causal-variant information that SuSiEx exploits.

### Case-control GWAS passing Ntotal instead of Neff

**Trigger:** Passing `n = N_total` to `susie_rss()` for case-control GWAS derived from logistic regression.

**Mechanism:** susie_rss expects the effective sample size that determined the standard errors. For case-control logistic regression, `Neff = 4 / (1/Ncase + 1/Ncontrol)`; when cases are rare, total N can exceed Neff by 25x or more. Passing Ntotal rescales z-scores into a regime SuSiE never sees and makes the implied prior variance wrong.

**Symptom:** PIPs systematically biased; credible sets either too narrow (PIPs collapse to a single SNP that is not robust) or too wide (PIPs flatten); replication poor; sometimes z-score scale warnings from susieR.

**Fix:** `Neff = 4 / (1/Ncase + 1/Ncontrol)`. Example: Ncase=5000, Ncontrol=495000 -> Neff ~= 19,800 (NOT 500,000). For quantitative traits from linear regression, `n = N_total` is correct. Reference: Privé 2022 HGG (`bigsnpr` documents Neff handling); Willer 2010 Bioinformatics (METAL Neff convention).

### Allele Harmonization with the LD Reference

**Trigger:** Effect allele in GWAS sumstats differs from the coding/A1 allele in the LD reference panel; or palindromic SNPs (A/T, C/G) carried without strand resolution.

**Mechanism:** susie_rss treats `z` and `R` as defined on the same allele coding. If the effect allele is swapped relative to the LD-reference A1, the sign of z is wrong and the LD row/column for that SNP is implicitly flipped. SNPs matching by rsID can silently swap alleles between sumstats and reference, breaking the `z' R z` consistency the model relies on.

**Symptom:** `estimate_s_rss` lambda inflated despite ancestry-matched panel; `kriging_rss` flags many SNPs with `|z_obs - z_exp| > 3` clustered at SNPs where reference A1 != GWAS effect allele; credible sets pick up tag-only SNPs anti-correlated with the lead.

**Fix:** Harmonize before fitting:

```r
harmonize_z_to_ref <- function(z, gwas_a1, gwas_a2, ref_a1, ref_a2) {
    palindromic <- (gwas_a1 == 'A' & gwas_a2 == 'T') | (gwas_a1 == 'T' & gwas_a2 == 'A') |
                   (gwas_a1 == 'C' & gwas_a2 == 'G') | (gwas_a1 == 'G' & gwas_a2 == 'C')
    flip <- (gwas_a1 == ref_a2) & (gwas_a2 == ref_a1)
    z[flip] <- -z[flip]
    drop <- palindromic | !((gwas_a1 == ref_a1 & gwas_a2 == ref_a2) | flip)
    list(z = z, keep = !drop)
}
```

Drop palindromic SNPs at MAF > 0.42 (ambiguous strand); or resolve via external strand info (TopMed, 1000G strand files). `TwoSampleMR::harmonise_data()` offers an alternative implementation. See causal-genomics/colocalization-analysis for an equivalent harmonize helper used downstream.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| SuSiE finds 3 credible sets, FINEMAP finds 1 | FINEMAP's stochastic search did not converge OR SuSiE absorbed background into spurious sets | Increase FINEMAP `--n-iter`; check SuSiE purity (sets with purity < 0.5 are spurious) |
| SuSiE PIPs much sharper than FINEMAP | susie_rss assumes single residual variance; FINEMAP marginalizes over noise | Both can be correct; report the intersection of high-PIP variants from both as primary candidates |
| PolyFun + SuSiE collapses 10-variant credible set to 1 | Functional priors are doing real work (coding variant in set) | Verify with `prior_weights` plot; if priors are coding-specific the result is interpretable |
| SuSiEx credible set excludes the EUR top-PIP variant | EUR signal is tag, true causal shared across ancestries lies elsewhere | Trust SuSiEx if both populations have well-powered GWAS; verify with conditional analysis |
| HLA gives 50-variant credible set in every method | HLA LD structure cannot be fine-mapped by linear methods | Use HLA-specific imputation (HIBAG, SNP2HLA) and haplotype-level analysis |

**Operational rule:** For high-confidence reporting, require that (a) `estimate_s_rss()` lambda < 0.05; (b) at least one credible set has purity > 0.5 (`min_abs_corr >= 0.5`, equivalent to r2 >= 0.25); (c) the lead PIP variant within that set is reproduced by an independent method (FINEMAP, SuSiEx, or in-sample SuSiE if reference-LD was used). Anything failing these three is exploratory.

## Quantitative Thresholds

| Quantity | Threshold | Source / Rationale |
|----------|-----------|-------------------|
| Credible set coverage (well-powered GWAS) | 0.95 (default) | Wang 2020 JRSSB; standard convention |
| Credible set coverage (low-N eQTL, GTEx tissue) | 0.9 or 0.8 | At N < 1000, default 0.95 returns very wide CS; document choice in methods |
| Credible set purity (rare-variant fine-mapping) | min_abs_corr >= 0.1 | LD genuinely sparse; relax to retain signal |
| Credible set purity (default common-variant) | min_abs_corr >= 0.5 (r2 >= 0.25) | susieR default; below this the set is LD-confounded |
| Credible set purity (publication-strict) | min_abs_corr >= 0.7 | Stringent claim; rare in practice |
| PIP suggestive | > 0.5 | Convention; "more likely than not causal among set" |
| PIP strong | > 0.9 | Convention; high-confidence single candidate |
| PIP very strong | > 0.95 | Convention; near-certain candidate within credible set |
| L (default cap) | 10 | susieR default; sufficient for most non-HLA loci |
| L (HLA / complex loci) | 20-30 | Empirical; HLA hosts > 10 independent signals for many traits |
| `n` for case-control susie_rss | Neff = 4/(1/Ncase + 1/Ncontrol), NOT Ntotal | Privé 2022 HGG; matches the SE scale of logistic-regression sumstats |
| `estimate_s_rss` lambda acceptable | < 0.05 | susieR vignette; > 0.10 indicates serious LD mismatch |
| `kriging_rss` per-SNP flag | |z_obs - z_exp| > 3 | susieR vignette; flag for manual review |
| Locus window (default) | +/- 500 kb from sentinel | Conventional; covers most LD blocks |
| Locus window (conditional-p floor) | Extend until conditional -log10(p) < 4 | Avoids truncating a secondary signal whose conditional evidence leaks into the window edge |
| Locus window (long-range LD) | 5+ Mb or stratify | HLA chr6:25-35Mb, chr8 inversion chr8:8.1-11.9Mb hg38, chr17 H1/H2 inversion |
| FINEMAP `--n-causal-snps` | 5 | Default; raise for HLA |
| FINEMAP `--sss-tol` | 0.001 | Convergence tolerance; rarely needs change |

## Functional Priors with PolyFun

**Goal:** Use genome-wide stratified LDSC heritability to weight per-SNP causal priors, sharpening PIPs at coding, conserved, and regulatory variants.

**Approach:** Run PolyFun once genome-wide to estimate per-SNP h2 from the baseline-LF annotation set; extract per-SNP causal prior; pass to susie_rss as `prior_weights`.

```bash
# Stage 1: compute per-SNP h2 contributions genome-wide (one-time)
polyfun.py --compute-h2-L2 \
    --output-prefix polyfun_h2 \
    --sumstats gwas_munged.sumstats \
    --ref-ld-chr UKB_baseline_LF/baselineLF2.2.UKB. \
    --w-ld-chr UKB_baseline_LF/weights.UKB.

# Stage 2: extract per-SNP prior for one locus
polyfun.py --compute-h2-bins \
    --output-prefix polyfun_priors \
    --sumstats gwas_munged.sumstats \
    --ref-ld-chr polyfun_h2
```

```r
library(susieR)
priors <- read.table('polyfun_priors.6.snpvar.gz', header = TRUE)
priors <- priors[match(gwas_df$SNP, priors$SNP), ]
prior_w <- priors$SNPVAR / sum(priors$SNPVAR, na.rm = TRUE)

fit <- susie_rss(z = z_scores, R = ld_matrix, n = N, L = 10,
                 prior_weights = prior_w)
```

UKB baseline-LF priors are pre-computed EUR-only at `data.broadinstitute.org/alkesgroup/UKBB_LD/` for hg19 and hg38. For EAS, AFR, or SAS GWAS, the EUR weights are NOT valid: use PolyFun-pop (Weissbrod 2022 Nat Genet 54:450) which extends PolyFun via cross-population annotation weighting, OR the cross-ancestry baseline-LF tables from Marquez-Luna 2025. For ancestries lacking matched baseline-LF (admixed, under-represented), accept reduced power and run uniform-prior susie_rss; applying EUR weights to non-EUR sumstats produces miscalibrated PIPs that look sharper than reality.

### Manual Coding-Variant Priors Without PolyFun

For postdocs without PolyFun infrastructure or with single-locus inputs, manual annotation-based priors are a reasonable approximation (Hutchinson 2020 Genet Epidemiol 44:313). Coding variants get ~10x uniform weight; broadly conserved variants ~5x (binned by CADD-PHRED quantile).

```r
build_manual_priors <- function(vep_df, cadd) {
    w <- rep(1, nrow(vep_df))
    w[vep_df$Consequence %in% c('missense_variant', 'stop_gained', 'splice_donor_variant',
                                'splice_acceptor_variant', 'frameshift_variant')] <- 10
    w[cadd >= quantile(cadd, 0.95, na.rm = TRUE)] <- pmax(w[cadd >= quantile(cadd, 0.95, na.rm = TRUE)], 5)
    w / sum(w)
}
fit <- susie_rss(z = z_scores, R = ld_matrix, n = Neff, L = 10, prior_weights = build_manual_priors(vep, cadd))
```

Report the prior construction explicitly; reviewers will ask whether the prior was tuned post hoc.

## Cross-Ancestry Fine-Mapping with SuSiEx

**Goal:** Jointly fine-map a locus across multiple ancestries assuming shared causal variants but population-specific LD.

**Approach:** Per-ancestry summary statistics + per-ancestry LD reference; SuSiEx runs a joint SuSiE model with population-specific R matrices.

```bash
SuSiEx \
    --sst_file=eur_sumstats.txt,eas_sumstats.txt,afr_sumstats.txt \
    --n_gwas=500000,200000,80000 \
    --ref_file=1000G_EUR,1000G_EAS,1000G_AFR \
    --ld_file=eur_ld,eas_ld,afr_ld \
    --out_dir=susiex_out \
    --out_name=locus1 \
    --pop=EUR,EAS,AFR \
    --chr=6 --bp=30000000,31000000 \
    --chr_col=1,1,1 --snp_col=2,2,2 --bp_col=3,3,3 \
    --a1_col=4,4,4 --a2_col=5,5,5 --eff_col=6,6,6 \
    --se_col=7,7,7 --pval_col=8,8,8 \
    --level=0.95
```

The output includes per-population PIPs and a joint credible set. Credible sets from SuSiEx are typically 2-5x smaller than EUR-only susie_rss when AFR is included, because AFR shorter LD blocks resolve EUR-tagged regions.

## FINEMAP CLI Pattern

**Goal:** Independent confirmation via shotgun stochastic search.

**Approach:** Build .z, .ld, and master files; run FINEMAP with `--sss` and parse the .snp and .cred outputs.

```bash
# .z file format: snp chromosome position allele1 allele2 maf beta se
# .ld file: square LD matrix, space-separated, no header

cat > locus.master <<'EOF'
z;ld;snp;config;cred;log;n_samples
locus.z;locus.ld;locus.snp;locus.config;locus.cred;locus.log;500000
EOF

finemap --sss \
    --in-files locus.master \
    --n-causal-snps 5 \
    --sss-tol 0.001 \
    --n-iter 100000 \
    --n-conv-sss 5000

# Parse:
# locus.snp -> per-variant prob (PIP), log10bf
# locus.cred -> credible sets at increasing causal counts
# locus.config -> top configurations
```

FINEMAP and SuSiE agree when sparsity holds; disagreement often reveals non-sparse loci that need SuSiE-inf.

## Coloc.susie Integration

**Goal:** Test colocalization between two traits using credible sets, not single SNPs.

**Approach:** Fit susie_rss separately per trait; pass both `susie` objects to `coloc.susie`; per-credible-set colocalization probabilities are returned.

```r
library(coloc)

fit_trait1 <- susie_rss(z = z1, R = ld_matrix, n = N1, L = 10)
fit_trait2 <- susie_rss(z = z2, R = ld_matrix, n = N2, L = 10)

coloc_res <- coloc.susie(fit_trait1, fit_trait2)
# coloc_res$summary: per-credible-set PP.H4 (shared causal probability)
print(coloc_res$summary)
```

PP.H4 > 0.8 per credible set is the conventional shared-causal threshold; weaker thresholds suggest distinct or conditional signals. See causal-genomics/colocalization-analysis.

## HLA and Long-Range LD: When to Stop

The HLA region (chr6:28-34 Mb), chromosome 8 inversion (chr8:8-12 Mb), and a handful of other extended LD blocks violate the assumptions of every fine-mapping method.

**Symptoms of irrecoverable LD structure:** Credible sets contain 30+ SNPs at low purity even with L=30; SuSiE-inf credible sets remain wide; `kriging_rss` flags hundreds of SNPs.

**Options:**
- Stratify by classical HLA allele (HIBAG, SNP2HLA imputation) and test allelic series
- Conditional analysis on the lead variant before fine-mapping the residual
- Exclude the region from genome-wide fine-mapping summaries and report separately
- For chr8 inversion: stratify by inversion genotype if known

Document the caveat in any methods section; standard PIPs at HLA are not interpretable as causality estimates.

## TWAS Fine-Mapping (FOCUS) -- delegated

Variant-level fine-mapping cannot distinguish causal genes among co-regulated TWAS hits. FOCUS / MA-FOCUS extend fine-mapping to the predicted-expression level; see causal-genomics/transcriptome-wide-association for the full FOCUS workflow and reconciliation with variant-level credible sets.

## Required Reporting Schema for Fine-Mapping

Every locus reported should carry these columns; missing fields are the most common reviewer complaint.

| Column | Description |
|--------|-------------|
| locus_id | Locus identifier (chr:start-end or sentinel rsID) |
| method | susie_rss / FINEMAP / PAINTOR / SuSiEx / SuSiE-inf |
| L_used | `sum(!fit$sets$pruned)` (effective L; not just the cap passed in) |
| n_credible_sets | Number of returned credible sets at the chosen coverage |
| cs_size | Variants per credible set |
| cs_purity_min / cs_purity_mean | min and mean `fit$sets$purity[,'min.abs.corr']` per set |
| top_pip_snp | Lead variant in each credible set |
| top_pip | Posterior inclusion probability of top_pip_snp |
| lambda_s | `estimate_s_rss` diagnostic for the locus |
| kriging_outlier_count | Count of SNPs with `|z_obs - z_exp| > 3` |
| ld_panel | 1KG-EUR / UKB-EUR / in-sample / TopMed |
| prior_source | uniform / PolyFun-EUR / PolyFun-pop / manual coding-variant |
| coverage | 0.95 default; 0.9 or 0.8 documented for low-N |
| n_effective | Sample size passed to susie_rss (Neff for case-control) |

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "In-sample vs reference LD?" | In-sample preferred when cohort genotypes available; if reference, report `estimate_s_rss` lambda < 0.05 plus `kriging_rss` outlier count |
| "Credible-set purity?" | `min_abs_corr >= 0.5` (r2 >= 0.25) default; reported per set; relaxed only with explicit rationale for rare-variant fine-mapping |
| "Is L set high enough?" | If returned CS count < L cap: OK (susieR auto-prunes); otherwise raise L. HLA needs L=20-30 |
| "Why not SuSiE-inf?" | Polygenic-shoulder test: count SNPs with marginal -log10(p) > 4 outside the lead credible set; > 50 indicates a polygenic shoulder and SuSiE-inf (Cui 2024) should be used |
| "Why no functional priors?" | PolyFun applied (or manual coding-variant prior used) and reported; if uniform, justify (low-N, mismatched-ancestry baseline-LF) |
| "Credible set has 50 SNPs -- is that fine-mapping?" | Acknowledged as imprecise; reported alongside diagnostics; cross-trait colocalization or functional fine-mapping (PolyFun, MPRA, allelic series) recommended for resolution |
| "Did you use Neff for case-control?" | Yes: `Neff = 4/(1/Ncase + 1/Ncontrol)`; report the value used |
| "Allele harmonization?" | Yes: flipped z when GWAS effect allele differs from reference A1; palindromic SNPs at MAF > 0.42 dropped |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `IBSS algorithm did not converge` warning | L too small OR LD mismatch | Increase L; run `estimate_s_rss`; check ancestry match |
| Credible set contains all SNPs at locus | LD reference matches discovery poorly; `s_hat` > 0.1 | Switch to in-sample LD or stratify reference ancestry |
| PIPs identical to GWAS p-value rank | Effectively no LD information used; check LD matrix orientation | Verify SNP order in z and R match exactly; check for transposed R |
| Negative eigenvalues in LD matrix | Numerical PSD violation from finite-precision storage | Add small ridge: `R <- R + diag(1e-4, nrow(R))`; or use `Matrix::nearPD` |
| `pip` all ~ 1/p (uniform) | Convergence failure OR all effects pruned | Check `fit$converged`; raise L; check Z scale |
| FINEMAP `Error: SNP names do not match` | .z and .ld SNP order differ | Ensure both are sorted identically; pass matched .snp file |
| Coloc.susie returns NULL | One trait has zero credible sets | Verify both fits succeeded; lower coverage to 0.9 if signal is weak |
| SuSiEx output empty | Population label mismatch with reference panel | Verify `--pop` matches reference filename convention; check `--bp` window |
| PolyFun priors do not change PIPs | Passed to `prior_variance` instead of `prior_weights` | Read susieR docs; use `prior_weights=` |

## References

- Wang G, Sarkar A, Carbonetto P, Stephens M 2020 J R Stat Soc B 82:1273 (SuSiE / IBSS)
- Zou Y, Carbonetto P, Wang G, Stephens M 2022 PLoS Genet 18:e1010299 (susie_rss for summary statistics)
- Cui R, Elzur RA, Kanai M, Ulirsch JC, Weissbrod O et al 2024 Nat Genet 56:162 (SuSiE-inf / FINEMAP-inf for non-sparse loci)
- Benner C, Spencer CC, Havulinna AS, Salomaa V, Ripatti S, Pirinen M 2016 Bioinformatics 32:1493 (FINEMAP)
- Hormozdiari F, Kostem E, Kang EY, Pasaniuc B, Eskin E 2014 Genetics 198:497 (CAVIAR)
- Wen X, Lee Y, Luca F, Pique-Regi R 2016 AJHG 98:1114 (DAP-G)
- Kichaev G, Yang WY, Lindstrom S, Hormozdiari F, Eskin E et al 2014 PLoS Genet 10:e1004722 (PAINTOR)
- Weissbrod O, Hormozdiari F, Benner C, Cui R, Ulirsch J et al 2020 Nat Genet 52:1355 (PolyFun + functional priors)
- Yuan K, Longchamps RJ, Pardinas AF, Yu M, Chen TT et al 2024 Nat Genet 56:1841 (SuSiEx cross-ancestry)
- Tashman KC, Cui R, O'Connor LJ, Pasaniuc B, Price AL 2024 medRxiv / 2025 Nat Genet (MultiSuSiE; published in Nat Genet 2025, doi:10.1038/s41588-025-02450-5)
- Mancuso N, Freund MK, Johnson R, Shi H, Kichaev G et al 2019 Nat Genet 51:675 (FOCUS for TWAS fine-mapping)
- Wallace C 2021 PLoS Genet 17:e1009440 (coloc.susie integration)
- Schaid DJ, Chen W, Larson NB 2018 Nat Rev Genet 19:491 (fine-mapping review)
- Hutchinson A, Asimit J, Wallace C 2020 Genet Epidemiol 44:313 (fine-mapping methods comparison)

## Related Skills

- causal-genomics/colocalization-analysis - coloc.susie consumes susie_rss credible sets; equivalent harmonize helper
- causal-genomics/effector-gene-prioritization - Downstream gene-assignment from credible-set variants
- causal-genomics/transcriptome-wide-association - FOCUS / MA-FOCUS for gene-level fine-mapping
- causal-genomics/genomic-sem - Joint multi-trait fine-mapping when credible sets are shared across traits
- causal-genomics/mendelian-randomization - Fine-mapped variants as cis-instruments
- causal-genomics/pleiotropy-detection - Per-credible-set pleiotropy testing
- atac-seq/enhancer-gene-linking - ABC / ENCODE-rE2G linking credible-set variants to target genes
- population-genetics/linkage-disequilibrium - Constructing LD matrices for susie_rss
- population-genetics/association-testing - Upstream GWAS summary statistic generation
- workflows/gwas-pipeline - End-to-end GWAS pipeline producing fine-mapping input
- variant-calling/variant-annotation - Annotating credible-set variants with VEP / coding consequence
- pathway-analysis/go-enrichment - Downstream gene-level interpretation of credible-set targets
