---
name: bio-causal-genomics-colocalization-analysis
description: Test whether two or more traits share a causal variant at a locus using Bayesian colocalization (coloc.abf, coloc.susie, HyPrColoc, moloc, eCAVIAR, SMR/HEIDI, PWCoCo, SharePro). Use when integrating GWAS with eQTL/sQTL/pQTL/mQTL, distinguishing shared causal variants from LD-driven coincidence, handling allelic heterogeneity, choosing between single-causal vs multi-causal methods, picking PP.H4 thresholds, running sensitivity over p12, or harmonising summary statistics for colocalization.
tool_type: r
primary_tool: coloc
---

## Version Compatibility

Reference examples tested with: coloc 5.2.3+, susieR 0.12.35+, hyprcoloc 1.0+ (GitHub jrs95/hyprcoloc), SMR 1.3.1+ (CLI, cnsgenomics.com), eCAVIAR 2.2+ (compiled from caviar/eCAVIAR repo), PWCoCo 1.0+ (jwr-git/pwcoco), moloc 0.1+ (clagiamba/moloc), SharePro_coloc 7.0+ (zhwm/SharePro_coloc), R >= 4.1.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('coloc')`; check `?coloc.abf`, `?coloc.susie`, `?runsusie`
- CLI: `smr --version`, `pwcoco --help`, `sharepro_coloc.py --help`

If code throws AttributeError, NULL list elements, or `Error in coloc.abf: dataset must have...`, introspect the installed package signature and adapt the example rather than retrying.

# Colocalization Analysis

**"Test whether my GWAS signal and an eQTL share the same causal variant"** -> Compute Bayesian posterior probabilities over five hypotheses (H0 neither, H1 trait-1-only, H2 trait-2-only, H3 distinct causal variants, H4 shared causal variant) to discriminate true causal overlap from LD-driven coincidence, then run sensitivity analysis over the p12 prior.

- R (single-causal, fastest): `coloc::coloc.abf(dataset1, dataset2, p12=5e-6)` -> `coloc::sensitivity(res, 'H4 > 0.75')`
- R (multi-causal, needs LD): `runsusie(d1)` -> `runsusie(d2)` -> `coloc.susie(s1, s2)` -> per-credible-set PP
- R (many traits, single-causal cluster): `hyprcoloc::hyprcoloc(effect.est = betas_mat, effect.se = ses_mat, trait.names = ..., snp.id = ...)` -> trait clusters
- CLI (causality vs linkage): `smr --bfile ref --gwas-summary g.ma --beqtl-summary eqtl.besd --out smr` -> SMR p + HEIDI p
- CLI (allelic heterogeneity): eCAVIAR `eCAVIAR -l ld1 -l ld2 -z z1 -z z2 -o out -c 2` -> CLPP per SNP
- CLI (conditional): PWCoCo conditions on each independent signal via GCTA-COJO then runs pairwise coloc.abf

## Algorithmic Taxonomy

| Method | Model | Inputs | Output | Strength | Fails when |
|--------|-------|--------|--------|----------|------------|
| coloc.abf (Giambartolomei 2014) | Single causal variant per locus; Bayesian ABF | beta+varbeta or p+MAF; sample sizes; type/s/sdY | PP.H0-H4 | Fast (~1s/locus), no LD required, mature, widely-cited | 2+ causal variants in moderate LD -> PP.H3 inflates spuriously; assumes a single causal per trait |
| coloc.susie (Wallace 2021) | Multi-causal via SuSiE; per-credible-set pairwise coloc | Summary stats + ancestry-matched LD matrix | PP.H4 per (CS1, CS2) pair | Handles allelic heterogeneity; principled CS framework | Sensitive to LD-mismatch; sample-size-LD mismatch -> spurious credible sets; needs in-sample or matched LD |
| SMR + HEIDI (Zhu 2016) | Tests pleiotropy (one variant -> both traits) vs linkage (two variants in LD) | GWAS .ma; eQTL .besd; LD reference (plink bfile) | SMR p (significance) + HEIDI p (null = shared causal) | Distinguishes shared-causal from linkage at a top SNP; standard for eQTLGen / GTEx integration | Fails to discriminate when LD between causal SNPs > 0.7 (HEIDI loses power); HEIDI requires >= 10 SNPs near top |
| eCAVIAR / CLPP (Hormozdiari 2016) | Fine-mapping-aware; computes Colocalization Posterior Probability per SNP | Z-scores; LD matrices per trait | CLPP per SNP; per-locus sum | Handles allelic heterogeneity natively; per-SNP resolution | Computationally heavy at -c > 3 causal variants; CLPP thresholds debated (0.01 vs 0.1) |
| PWCoCo (Robinson 2022) | Pairwise conditional via GCTA-COJO conditioning | Summary stats + individual-level LD bfile | Per-conditional-signal coloc.abf results | Cleanly handles AH at top GWAS hit + secondary signals | Needs individual-level reference; sensitive to COJO collinearity threshold |
| moloc (Giambartolomei 2018) | Multi-trait extension of coloc.abf (3-5 traits) | Per-trait summary stats | 15 (3-trait) / 31 (4-trait) / 63 (5-trait) hypothesis PPs | First principled multi-omic coloc | Hypothesis count = 2^k - 1 explodes; >= 6 traits computationally infeasible; minimally updated since 2019 |
| HyPrColoc (Foley 2021) | Many-trait cluster-based; iterative branch-and-bound under single-causal | Beta + SE matrices SNPs x traits | Trait clusters sharing a causal variant | Scales to 50+ traits; identifies cluster substructure | Inherits single-causal assumption from coloc.abf; clusters can fragment under AH |
| SharePro_coloc (Wang 2024) | Variational effect-group joint model | Beta + SE; LD per ancestry | Effect-group level PP | Handles multi-causal + multi-ancestry jointly; faster than coloc.susie at scale | Newer (2024); benchmarks evolving; trickier installation |
| Wallace 2020 / 2025 variant-specific priors | Function-aware p12 (e.g. up-weight coding/promoter SNPs) | Same as coloc.abf + per-SNP prior weights | PP.H0-H4 with non-uniform prior | Improves discovery when functional annotation is informative | Annotation choice is a methodological lever; report sensitivity |

Methodology evolves; verify the current Open Targets Genetics, eQTL Catalogue, and FinnGen colocalization pipelines before locking parameters. Open Targets uses coloc.abf at PP.H4 >= 0.75 with p12 = 1e-5; FinnGen uses coloc.susie at PP.H4 >= 0.8 with in-sample LD.

## Decision Tree by Scenario

| Scenario | Recommended method | Why |
|----------|---------------------|-----|
| GWAS + single-tissue eQTL, top GWAS variant looks single-signal | coloc.abf + sensitivity() | Fast, no LD needed, well-validated; single-causal assumption typically holds at clean loci |
| GWAS + eQTL, conditional analysis shows 2+ independent signals | coloc.susie OR PWCoCo | Multi-causal handling; coloc.susie if summary-stats LD available, PWCoCo if individual-level reference accessible |
| GWAS + multi-tissue eQTL (e.g. all 49 GTEx tissues) | coloc.abf per tissue + HyPrColoc across tissues | Per-tissue PP.H4 gives tissue-specific causality; HyPrColoc identifies tissue clusters sharing the variant |
| GWAS + eQTL + sQTL + mQTL (3-5 omics) | moloc (k <= 5) OR HyPrColoc | moloc gives explicit hypothesis posterior; HyPrColoc scales but loses hypothesis structure |
| Many GWAS traits at one locus (pleiotropic hub) | HyPrColoc | Designed for many-trait clustering; coloc.abf pairwise scales as k^2 |
| Top SNP has only modest GWAS p; is it the same causal as eQTL? | SMR + HEIDI | SMR tests pleiotropy/causality; HEIDI rejects shared-causal -> linkage |
| Want per-SNP credibility under allelic heterogeneity | eCAVIAR (CLPP) | Per-SNP CLPP integrates fine-mapping with coloc |
| MHC / HLA region (chr6:25-35 Mb) | HLA-coloc (Lagou 2024) OR exclude MHC | Long-range LD breaks single-causal assumption; standard PP.H4 not interpretable |
| Trans-eQTL / GWAS pair | coloc.abf with p12 lowered to 5e-6 or 1e-6 | Shared causality is biologically rare; default p12=1e-5 over-favours H4 |
| Ancestry-mismatched GWAS vs eQTL | SharePro_coloc OR ancestry-matched coloc.susie | LD differs across ancestries; using EUR LD on AFR z-scores produces spurious credible sets |
| Very small eQTL (N < 200) | None reliably; flag locus underpowered | All methods report H0/H1/H2 dominance; report PP transparently and gather larger reference (eQTLGen N~31k, GTEx v8) |

## Per-Method Failure Modes

### coloc.abf -- PP.H3 inflation under multiple causal variants

**Trigger:** Locus has 2+ independent causal signals in moderate LD (r2 ~ 0.3-0.6).

**Mechanism:** The single-causal-variant assumption forces the model to allocate posterior mass to H3 (distinct causal variants) whenever the per-SNP Bayes factors for the two top SNPs do not align.

**Symptom:** Visual co-localization in LocusZoom looks convincing, but `result$summary['PP.H3.abf']` dominates over PP.H4; sensitivity() shows PP.H4 stays low across the entire p12 grid.

**Fix:** Run coloc.susie (or eCAVIAR or PWCoCo) to allow multiple causal variants. If coloc.susie returns multiple credible sets with one pair showing PP.H4 > 0.75, this is real allelic heterogeneity not failure.

### coloc.susie -- LD reference mismatch

**Trigger:** Z-scores from GWAS / eQTL of one ancestry, LD matrix from 1000 Genomes EUR (or any non-matched reference).

**Mechanism:** SuSiE assumes z-scores and the supplied LD are jointly consistent. Ancestry mismatch or sample-size mismatch produces a non-positive-definite implicit covariance; SuSiE responds by returning spurious credible sets that include LD-mismatched SNPs.

**Symptom:** `susieR::estimate_s_rss(z, R, n)` returns lambda > 0.05; `susieR::kriging_rss` flags off-diagonal SNPs with extreme studentized residuals; credible sets are oddly large (50+ SNPs) or include SNPs distant in LD from the lead.

**Fix:** Use in-sample LD when at all possible (per-cohort plink `--r square`). If reference must be external, match ancestry (1KG superpopulation) and superpopulation-stratify. Run `estimate_s_rss` and report lambda; if > 0.05, drop the locus or switch to coloc.abf.

### coloc default p12 too liberal for trans-eQTL

**Trigger:** Applying `p12 = 1e-5` (the default) to a trans-eQTL / GWAS pair.

**Mechanism:** The default p12 was calibrated for cis-eQTL where biological proximity makes shared causality reasonable. For trans-eQTL, prior probability of shared causality is much lower; uniform p12 over-favours H4.

**Symptom:** PP.H4 > 0.8 reported, but sensitivity() reveals PP.H4 falls below 0.5 for p12 < 1e-5; replication in independent data fails.

**Fix:** Operational definition: "trans" = >5 Mb from TSS or different chromosome. Default p12=1e-5 over-favors H4 for trans (genome-rare biology). For trans: lower p12 to 5e-6 or 1e-6 AND raise PP.H4 threshold to >= 0.8 (compensate for higher FP risk). Cross-reference Vosa 2021 Nat Genet 53:1300 (eQTLGen trans) for empirical patterns.

### MHC / HLA + chr 8 inversion -- single-causal assumption breaks

**Trigger:** Locus within chr6:25-35 Mb (extended MHC, hg38), or chr8:8.1-11.9 Mb (chr 8 inversion, hg38).

**Mechanism:** The MHC contains classical HLA genes with extreme long-range LD (r2 > 0.5 over many Mb), multiple independent causal haplotypes, and structural variation. The chr 8p23.1 inversion similarly produces long-range LD across megabases of polymorphic inversion alleles. The single-causal-variant assumption is biologically wrong in both regions.

**Symptom:** coloc.abf almost always returns PP.H3 or fragmented PP across H1/H2/H3/H4 even when the underlying biology is well-established (e.g. HLA-DRB1 in autoimmune GWAS).

**Fix (MHC):** Use HLA-imputed classical alleles via SNP2HLA / HIBAG / HLA-TAPAS, then HLA-coloc (Lagou 2024 medRxiv) -- NOT coloc on SNPs in MHC. OR exclude MHC from genome-wide coloc and report HLA association at the haplotype/allele level. **Fix (chr 8 inversion):** Exclude chr8:8.1-11.9 Mb or pre-condition on inversion genotype before coloc. Never report a single coloc PP.H4 in either region without this caveat.

### Lead-SNP swap and window bias

**Trigger:** The two traits have different lead SNPs at the same locus; analyst centres each window on the trait-specific lead.

**Mechanism:** coloc PP is sensitive to the SNPs in the window; centring on different leads gives different per-SNP overlap and biases toward H3.

**Symptom:** Re-centring the window on the GWAS lead vs the eQTL lead produces qualitatively different PP.H4.

**Fix:** Use a SINGLE window (typically +/- 500 kb or 1 Mb) centred on the joint top-variant (the SNP with the lowest min-p across both traits), or on the GWAS lead consistently. Report PP under multiple centring choices; flag the locus if PP swings > 0.2 across centrings.

### Underpowered eQTL (N < 200)

**Trigger:** Small eQTL discovery (e.g. tissue-specific bulk study, N < 200; per-cell-type sc-eQTL).

**Mechanism:** With low N, varbeta is large; the eQTL's per-SNP Bayes factors are flat; the joint likelihood concentrates on H0 or H1 (GWAS-only).

**Symptom:** PP.H0 or PP.H1 dominates; the eQTL panel shows visible signal but coloc cannot resolve causal vs noise.

**Fix:** Use eQTLGen (N ~ 31k whole-blood) or GTEx v8 (N ~ 70-700 per tissue) where possible. For rare cell types, accept the limitation and report the locus as underpowered rather than claim absence of colocalization.

| eQTL N | Coloc viability | Notes |
|--------|-----------------|-------|
| < 200 | Underpowered | PP.H1 dominant; flag |
| 200-500 | Cis only, modest | Single-tissue cis |
| 500-1000 | Good for cis | Most GTEx v8 tissues |
| >= 1000 | Well-powered | Trans accessible |
| >= 10000 | Meta (eQTLGen) | Cross-tissue / sc |

### Reference QTL panel choice

GTEx v8 (838 donors, 49 tissues, 2020) is the current PredictDB-supported standard. GTEx v10 (released 2024) has limited harmonisation and is not yet PredictDB-default. eQTLGen blood meta-eQTL (N ~ 31k) wins on sample size for blood cis-eQTL discovery, beating any single tissue on power. Always pin version in methods (e.g. "GTEx v8 MASHR-EUR, PredictDB release 2022-01").

## PP.H4 Threshold Framework

| Threshold | Use case | Source |
|-----------|----------|--------|
| 0.5 - 0.7 | Suggestive / pilot / hypothesis-generating | Giambartolomei 2014 original |
| >= 0.7 | Triangulation tier for TWAS / cis-MR / effector-gene cross-evidence | Open Targets Genetics common practice; cross-reference downstream skills |
| >= 0.75 | Open Targets Platform / eQTL Catalogue / FinnGen default screening threshold | Open Targets Genetics docs; Mountjoy 2021 Nat Genet 53:1527 |
| >= 0.80 | Most published colocalizations / standard publication tier | Wallace 2020 PLoS Genet 16:e1008720 |
| >= 0.90 | Stringent clinical / therapeutic-target prioritization | Reserved for high-confidence claims |
| >= 0.95 | Industry / regulatory drug-target submission grade | Internal pharma default |
| PP.H3 >= 0.80 | Confident distinct causal variants (negative coloc result) | Standard |
| PP.H4 / (PP.H3 + PP.H4) >= 0.9 | Conditional probability framing (some pipelines) | Foley 2021 |

**Operational rule:** Three operational tiers map onto the most common downstream uses: (a) **>= 0.7** when PP.H4 is one of several lines of triangulating evidence (TWAS + coloc, cis-MR + coloc, effector-gene multi-evidence) -- this is the threshold downstream skills (causal-genomics/transcriptome-wide-association, causal-genomics/mendelian-randomization cis-MR, causal-genomics/effector-gene-prioritization, causal-genomics/proteome-mr-drug-target) require; (b) **>= 0.8** for standard peer-reviewed publication as a stand-alone coloc claim (Wallace 2020); (c) **>= 0.95** for industry / clinical drug-target submission. Open Targets and FinnGen pipelines screen at >= 0.75 but downstream publication-grade coloc claims should clear >= 0.8 and triangulation claims >= 0.7. ALWAYS report PP.H3 alongside PP.H4 -- a locus with PP.H4 = 0.6, PP.H3 = 0.3 is qualitatively different from PP.H4 = 0.6, PP.H3 = 0.05 (the former is real ambiguity over single vs distinct causal; the latter is underpowered evidence). Run `coloc::sensitivity()` and report the p12 range over which PP.H4 stays above the threshold.

## Default Priors and the p12 Sensitivity Question

| Prior | Default | Interpretation | When to change |
|-------|---------|----------------|----------------|
| p1 | 1e-4 | Prob a random SNP is associated with trait 1 | Rarely changed |
| p2 | 1e-4 | Prob a random SNP is associated with trait 2 | Rarely changed |
| p12 | 1e-5 | Prob a random SNP is associated with both traits | Lower (5e-6 or 1e-6) for trans-eQTL or unrelated trait pairs; raise (5e-5) only with strong prior, e.g. molecular QTL in the same tissue as causal cell type |

The p12/p1 ratio (= 0.1 under defaults) is the prior odds of colocalization given a trait-1 association. Wallace 2020 (PLoS Genet 16:e1008720) showed default p12 = 1e-5 is too liberal for many real-world settings and recommended sensitivity analysis as standard practice. Wallace 2025 (PLoS Genet 21:e1011697) extended this with variant-specific priors weighted by functional annotation.

### p12 Sensitivity Grid

| p12 grid point | Use case | Reporting rule |
|----------------|----------|-----------------|
| 1e-4 | Suggestive only / EUR cis-eQTL relaxed | PP.H4 here cannot support a publication claim |
| 1e-5 | Default for most cis-eQTL <-> GWAS pairs | Standard |
| 5e-6 | Conservative cis; default for trans-eQTL coloc | Recommended publication baseline |
| 1e-6 | Very conservative; trans coloc with weak prior | Required for cross-trait genome-rare coloc |

**Operational rule:** Require PP.H4 to remain above threshold across at least 3 adjacent grid points; report the lowest p12 at which PP.H4 >= 0.75. Use `coloc::sensitivity(result, rule = 'H4 > 0.75')` for the diagnostic plot.

Required reporting: PP.H4 at default priors + p12 range over which PP.H4 stays above threshold.

## eCAVIAR CLPP Threshold Framework

CLPP (Colocalization Posterior Probability) is the per-SNP product of the two per-trait fine-mapping posteriors. Threshold conventions:

- Hormozdiari 2016 AJHG 99:1245 used CLPP >= 0.01 (validated against null simulations).
- 2024 GTEx / Open Targets pipelines use CLPP >= 0.05.
- High-confidence claims require CLPP >= 0.1.
- Report both sum-CLPP across the credible set AND max-CLPP at any single SNP -- the two answer different questions (locus-level vs lead-SNP-level confidence).
- LaPierre 2021 Bioinformatics (CAVIARBF) extends CLPP to conditional analysis.

```bash
eCAVIAR -l ld_gwas.ld -l ld_eqtl.ld \
        -z gwas.z -z eqtl.z \
        -o coloc_out -c 2     # -c = max independent causal variants per trait
# Output: per-SNP CLPP in coloc_out_col file; report sum and max
```

## LD Matrix Construction for coloc.susie

**Requirements:**

- Signed Pearson r (not r2). coloc.susie expects directional LD; squared LD silently inverts effect-direction inference.
- Ancestry-matched to GWAS / eQTL ancestry. EUR LD on AFR z-scores produces spurious credible sets.
- SNP-order-aligned to the beta vector and named to match (row/column names = SNP IDs).
- Positive semi-definite. Numerical-noise negative eigenvalues must be repaired.
- Effective N sample-size-matched to the trait being fine-mapped (provide via `runsusie(..., n = N)`).

```bash
# plink2 phased r (signed Pearson); square matrix output
plink2 --pfile 1KG_EUR \
    --extract snps.txt --chr 6 --from-bp X --to-bp Y \
    --r-phased square spaces --out locus_ld
```

```r
# Alternative: in-sample LD from BED via bigsnpr
R <- bigsnpr::snp_cor(snp_obj$genotypes, ind.col = locus_snps)
# PSD repair if negative eigenvalues from numerical noise
R <- as.matrix(Matrix::nearPD(R)$mat)
dimnames(R) <- list(snp_ids, snp_ids)
```

**Critical:** Row and column order of R MUST match SNP order in the beta vector -- silent failure otherwise. The SuSiE objective stays finite under mis-ordering and returns nonsense credible sets. Verify with `stopifnot(rownames(R) == names(beta))` before `runsusie`. Cross-reference causal-genomics/fine-mapping for the full LD diagnostic block (`estimate_s_rss` lambda < 0.05, `kriging_rss` outlier inspection).

## SMR vs coloc Reconciliation

SMR (Zhu 2016) and coloc test related but non-identical questions:

- **SMR** tests pleiotropy vs linkage: does the top eQTL SNP show a GWAS effect explainable by its eQTL effect (pleiotropic / causal) or does the GWAS effect come from a different SNP in LD (linkage)?
- **coloc** tests shared vs distinct causal variants over an entire window of SNPs.
- **HEIDI** is SMR's heterogeneity test; null hypothesis is single shared causal SNP. Zhu 2016 Nat Genet 48:481 specifies **HEIDI p > 0.05** (NOT 0.01) as non-rejection of single shared causal. HEIDI p > 0.05 does NOT prove shared causality -- only that data cannot reject it; pair with SMR p Bonferroni-corrected across probes. When LD between causal SNPs > 0.7, HEIDI loses power same as coloc.

When LD between two true causal SNPs is high (r2 > 0.7), both SMR/HEIDI and coloc.abf lose discriminatory power: SMR cannot pick which of the LD-tied SNPs is causal, and coloc.abf cannot reject H4 even if biology is two-distinct-causal. coloc.susie + ancestry-matched LD is the modern resolution.

**Operational rule:** SMR + HEIDI is appropriate when the question is "does this eQTL gene mediate the GWAS effect at all?" coloc is appropriate when the question is "do the two traits share a causal variant in this window?" Run both; agreement (significant SMR + non-rejected HEIDI + PP.H4 >= 0.75) is high-confidence; disagreement requires inspection (often the multi-causal / LD scenario above).

## moloc Multi-Omic Framework (3-5 Traits)

For k traits, moloc tests `2^k - 1` hypotheses. 3 traits -> 15 hypotheses (H_a, H_b, H_c, H_ab, H_ac, H_bc, H_abc, plus "none of the above"); 4 traits -> 31; 5 traits -> 63. The hypothesis H_{all-share} (all k share a single causal variant) is the multi-omic analog of PP.H4.

```r
# moloc 3-trait example; install via remotes::install_github('clagiamba/moloc')
library(moloc)
# Input: list of k dataframes with PVAL/N/F (frequency) per SNP and shared SNP IDs
result_moloc <- moloc_test(listData=list(gwas=gwas_df, eqtl=eqtl_df, sqtl=sqtl_df),
                            prior_var=c(0.01, 0.1, 0.5), priors=c(1e-4, 1e-6, 1e-7))
# PPA: posterior over all 15 hypotheses (3-trait case)
# Key column: PPA.abc (all-three-share)
```

moloc is computationally tractable up to k = 5 but explodes beyond; use HyPrColoc for k >= 6.

## HyPrColoc Cluster-Based Coloc (Many Traits)

HyPrColoc (Foley 2021) extends single-causal coloc to many traits by clustering traits that share a causal variant. Output: per-cluster posterior + per-trait cluster assignment.

```r
library(hyprcoloc)
# Inputs: SNPs-by-traits matrices of betas and standard errors
# Rows = SNPs (must be shared across all traits); Columns = traits
res <- hyprcoloc(effect.est=betas, effect.se=ses,
                  trait.names=colnames(betas), snp.id=rownames(betas),
                  reg.thresh=0.7,     # regional probability of coloc threshold
                  align.thresh=0.7)   # alignment threshold for traits within a cluster
res$results  # cluster assignment per trait + posterior
```

HyPrColoc inherits the single-causal-per-cluster assumption from coloc.abf; clusters can fragment if the true biology is allelic heterogeneity.

## PWCoCo (Conditional Pairwise Coloc)

PWCoCo (Robinson 2022) wraps GCTA-COJO conditional analysis around coloc.abf. For a locus with `k1` independent trait-1 signals and `k2` independent trait-2 signals, PWCoCo runs `k1 * k2` pairwise coloc.abf tests after conditioning each summary statistic on the other independent signals.

**When to use:** When GCTA-COJO has identified >= 2 independent signals in at least one trait and individual-level reference genotypes are available. Particularly suited to bulk eQTL with secondary cis signals.

**Inputs:** Per-trait summary stats (SNP, A1, A2, freq, beta, se, p, N) + plink bfile reference. **Output:** One coloc.abf result per (conditional signal 1, conditional signal 2) pair. Interpret each row as an independent single-signal coloc.

**Caveats:** PWCoCo requires individual-level reference (plink bfile); cannot run on summary stats alone. Collinearity threshold in COJO (default `--cojo-collinear 0.9`) controls how aggressively independent signals are split; lower values fragment, higher values merge. Worked CLI recipe in usage-guide.md.

## Standard coloc.abf Pipeline

**Goal:** Test whether a single GWAS lead variant shares a causal variant with an eQTL gene's top signal at a defined locus.

**Approach:** Extract a 1 Mb window centred on the GWAS lead; harmonise alleles between datasets; format coloc input lists with `type` ('quant' or 'cc'), sample size `N`, and either `sdY` (quant) or `s` (cc); run coloc.abf; run sensitivity() over the p12 grid.

```r
library(coloc)

# Inputs: harmonised gwas_df and eqtl_df with SNP, BETA, SE, MAF, N, POS columns
# Both must share the same SNP set and allele coding (verify with harmonise step)

gwas_input <- list(
    beta=gwas_df$BETA, varbeta=gwas_df$SE^2,
    snp=gwas_df$SNP, position=gwas_df$POS,
    type='cc',           # case-control GWAS
    s=0.30,              # case fraction
    N=50000)

eqtl_input <- list(
    beta=eqtl_df$BETA, varbeta=eqtl_df$SE^2,
    snp=eqtl_df$SNP, position=eqtl_df$POS,
    type='quant',        # quantitative eQTL
    sdY=1,               # SD(expression); 1 if standardised, else estimate from MAF+varbeta
    N=500)

res <- coloc.abf(dataset1=gwas_input, dataset2=eqtl_input,
                  p1=1e-4, p2=1e-4, p12=5e-6)   # conservative p12

print(res$summary)
sens <- coloc::sensitivity(res, rule='H4 > 0.75')   # generates plot + table
```

`sdY` semantics: when omitted, coloc estimates from `MAF` and `varbeta`; supplying `sdY=1` ASSUMES the trait is already standardised (eQTL with inverse-normal-transformed expression). Mismatch produces silently wrong Bayes factors -- the most common silent failure.

- For quantitative trait: leave `sdY=NULL` to estimate via `coloc:::sdY.est(varbeta, MAF, N)`. If CV of estimated sdY across SNPs > 0.5, varbeta/MAF are inconsistent -- coloc will silently miscalibrate Bayes factors.
- For inverse-normal-transformed expression: use `sdY = 1` (already standardized).
- Mismatch produces silently wrong PP -- most common silent failure.

`s` parameter for case-control (`type='cc'`):

- `s` = N_cases / N_total (NOT cases-per-control; NOT 0.5 default).
- For population-cohort case-control: `s` ~ disease prevalence in the cohort (~0.005 for rare disease).
- Wrong `s` does not error -- silently biases PP at extreme MAF.

## coloc.susie Multi-Causal Pipeline

**Goal:** Test colocalization at a locus with multiple independent signals (allelic heterogeneity).

**Approach:** Run SuSiE on each trait's summary statistics with ancestry-matched LD; verify LD-z-score consistency; coloc-test each pair of credible sets.

```r
library(coloc); library(susieR)

# Diagnostic: z-score vs LD consistency MUST be checked
z_gwas <- gwas_df$BETA / gwas_df$SE
lam_gwas <- susieR::estimate_s_rss(z=z_gwas, R=ld_matrix, n=gwas_n)
if (lam_gwas > 0.05) stop('LD reference mismatched to z-scores; lambda=', lam_gwas)

s1 <- runsusie(list(beta=gwas_df$BETA, varbeta=gwas_df$SE^2,
                    snp=gwas_df$SNP, position=gwas_df$POS,
                    type='cc', s=0.3, N=50000, LD=ld_matrix), L=10)
s2 <- runsusie(list(beta=eqtl_df$BETA, varbeta=eqtl_df$SE^2,
                    snp=eqtl_df$SNP, position=eqtl_df$POS,
                    type='quant', sdY=1, N=500, LD=ld_matrix), L=10)

res_susie <- coloc.susie(s1, s2)   # NULL if no overlapping CS
# res_susie$summary rows: each (hit1, hit2) pair of credible sets
```

LD matrix MUST be in the same SNP order as the beta vector; mis-ordering silently produces nonsense.

## SMR + HEIDI Pipeline

```bash
# SMR is a command-line tool. Pre-format GWAS into .ma (SNP A1 A2 freq beta se p N).
# eQTL data as BESD (binary eQTL summary data); pre-built BESD available from eQTLGen / GTEx.

smr --bfile 1KG_EUR_chr6 \
    --gwas-summary gwas.ma \
    --beqtl-summary eqtl_chr6.besd \
    --out smr_result \
    --thread-num 4 \
    --peqtl-smr 5e-8 \
    --heidi-mtd 1
# Output smr_result.smr: probe (gene) | top SNP | p_SMR | p_HEIDI | nsnp_HEIDI
```

Interpretation: significant `p_SMR` (Bonferroni-corrected across probes tested, typically < 5e-8 / N_probes) AND non-rejection by HEIDI (`p_HEIDI > 0.05`, per Zhu 2016) indicates pleiotropy / shared causal; `p_HEIDI <= 0.05` rejects shared-causal -> linkage. HEIDI p > 0.05 does NOT prove shared causality, only that data cannot reject it. Require `nsnp_HEIDI >= 10` for HEIDI reliability.

## Allele Harmonisation (Critical Pre-Step)

Mismatched effect alleles silently invert signs of betas, collapsing PP.H4 into PP.H3. Required steps before coloc:

1. Merge GWAS and eQTL summary stats by SNP ID (rsID or chr:pos:ref:alt).
2. Mark SNP-pairs as `same` (A1/A2 match) or `flip` (A1/A2 swap); drop SNPs that match neither.
3. For `flip` rows, negate the second dataset's beta (and swap A1/A2).
4. Drop palindromic SNPs (A/T or C/G) at MAF > 0.42; their strand cannot be inferred from coding alone (TwoSampleMR `harmonise_data` standard cutoff).
5. Verify genome build alignment (hg19 vs hg38 must match; lift over if not).

Worked harmonisation code and build-mismatch pitfalls: see usage-guide.md.

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Sensitivity to p12 prior?" | `coloc::sensitivity()` reported; PP.H4 robust across 1e-7 to 1e-5 grid |
| "Why not coloc.susie? Multi-causal possible?" | coloc.abf single-causal assumption stated; if PP.H3 dominant or GCTA-COJO identifies >= 2 independent signals, coloc.susie / SuSiE-based run; reported |
| "LD reference matched?" | In-sample preferred; if reference panel used, `estimate_s_rss(z, R, N)` lambda < 0.05; `kriging_rss` diagnostic clean |
| "PP.H4 = 0.6 is colocalization?" | No -- bands stated: 0.5-0.7 suggestive; >= 0.7 triangulation tier; >= 0.8 standard publication; >= 0.95 industry/clinical |
| "MHC region included?" | chr6:25-35 Mb excluded; HLA-coloc (Lagou 2024) for classical-allele-level coloc |
| "Ancestry mismatch?" | LD reference ancestry-matched to GWAS; for cross-ancestry use SharePro_coloc or coloc.susiex |
| "Sentinel SNP swap?" | Re-centered window on each trait's lead, joint top, eQTL top; PP.H4 stable within 0.1 |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `Error in coloc.abf: dataset must have N` | Forgot `N` in list, or `type` not set | Supply both; case-control also needs `s`; quant also needs `sdY` |
| PP.H3 dominant despite obvious visual overlap | 2+ causal in moderate LD breaks single-causal assumption | Run coloc.susie or eCAVIAR |
| `estimate_s_rss` lambda > 0.05 | LD reference does not match z-scores (ancestry / sample / build) | Use in-sample LD or ancestry-matched reference; do not proceed |
| PP.H4 unstable across p12 sensitivity grid | Borderline evidence; default priors not justified | Report the p12 range; lower priors for trans-eQTL; do not over-claim |
| coloc.susie returns NULL summary | No overlapping credible sets between traits | Genuine result (no shared signal) or both traits underpowered |
| Per-SNP betas have opposite signs but same magnitude across traits | Effect-allele mismatch | Run harmonisation; flip betas where A1/A2 swap; drop palindromic at high MAF |
| SMR significant + HEIDI p <= 0.05 | Linkage, not shared causal | Report as linkage; do not call colocalization |
| moloc all-share PPA collapses to ~0 | Different sample sizes / power across omics | Inspect per-omic effect sizes; consider HyPrColoc for cluster output |
| HyPrColoc trait cluster fragments | Underlying biology is multi-causal | Switch to coloc.susie at each suspected cluster centre |
| MHC PP.H4 close to 0 with strong visual signal | Long-range LD breaks single-causal | Use HLA-coloc or exclude MHC; never report standard coloc PP at MHC |

## Tool Install Notes

- **coloc**: CRAN. `install.packages('coloc')`. Bundles susieR dependency for >= 5.1.
- **susieR**: CRAN. `install.packages('susieR')`. >= 0.12.35 for `estimate_s_rss` and `kriging_rss`.
- **HyPrColoc**: GitHub only (never CRAN). `remotes::install_github('jrs95/hyprcoloc')`. Requires R >= 3.5.
- **SMR**: Pre-compiled binary from cnsgenomics.com/software/smr. Linux/Mac/Windows binaries; no R package.
- **eCAVIAR**: Compile from GitHub fhormoz/caviar; C++ source. CLI `eCAVIAR`. PAINTOR is the related multi-trait fine-mapping toolkit.
- **PWCoCo**: GitHub jwr-git/pwcoco. Compiled C++ CLI; can also be invoked from R via wrapper scripts.
- **SharePro_coloc**: Python `pip install sharepro-coloc` or GitHub zhwm/SharePro_coloc.
- **moloc**: GitHub clagiamba/moloc. R package; minimally updated since 2019, no CRAN release. R >= 3.5.

## Reviewer-Grade Reporting Template

For each colocalization claim, the report should include:

1. **Method and version** (e.g. coloc 5.2.3 coloc.abf, or coloc.susie with SuSiE L=10).
2. **Window definition** (e.g. +/- 500 kb around the GWAS lead rs12345 at chr6:30450000, hg38), and lead-SNP-swap sensitivity (PP.H4 at GWAS lead vs eQTL lead vs joint top).
3. **Priors** p1, p2, p12 used; **sensitivity** plot from `coloc::sensitivity()` and the p12 range over which PP.H4 stays above threshold.
4. **All five posteriors** PP.H0 through PP.H4 (not PP.H4 alone).
5. **Threshold band** the result clears (>= 0.7 triangulation tier, >= 0.75 Open Targets screening, >= 0.80 published, >= 0.90 stringent, >= 0.95 clinical).
6. **LD reference** ancestry, source (1000G phase 3 EUR / in-sample / UKBB), and lambda from `estimate_s_rss` if coloc.susie.
7. **Reference QTL panel** version (e.g. GTEx v8 MASHR-EUR, PredictDB release 2022-01; eQTLGen 2019).
8. **Conditional analysis** GCTA-COJO results if multi-causal; per-credible-set PP if coloc.susie.
9. **Failure-mode caveats** explicitly addressed: MHC excluded, chr 8 inversion excluded, ancestry-matched LD, sdY/s correctly specified, palindromic SNPs handled.
10. **Methods-section H0-H4 prose** describing what each hypothesis means (see usage-guide.md).

## References

- Giambartolomei C et al 2014 PLoS Genet 10:e1004383 (coloc.abf)
- Wallace C 2020 PLoS Genet 16:e1008720 (default-prior sensitivity; variant-specific priors)
- Wallace C 2025 PLoS Genet 21:e1011697 (further p12 prior refinements)
- Wallace C 2021 PLoS Genet 17:e1009440 (coloc.susie; multiple causal variants)
- Zhu Z et al 2016 Nat Genet 48:481 (SMR + HEIDI)
- Hormozdiari F et al 2016 AJHG 99:1245 (eCAVIAR / CLPP)
- Giambartolomei C et al 2018 Bioinformatics 34:2538 (moloc)
- Foley CN et al 2021 Nat Commun 12:764 (HyPrColoc)
- Robinson JW et al 2022 bioRxiv 2022.08.08.503158 (PWCoCo)
- Wang W et al 2024 Bioinformatics 40:btae295 (SharePro_coloc)
- Lagou V et al 2024 medRxiv (HLA-coloc)
- Mountjoy E et al 2021 Nat Genet 53:1527 (Open Targets Genetics colocalization pipeline)
- Vosa U et al 2021 Nat Genet 53:1300 (eQTLGen, N ~ 31,684 whole blood)
- GTEx Consortium 2020 Science 369:1318 (GTEx v8 multi-tissue eQTL)

## Related Skills

- causal-genomics/mendelian-randomization - Causal effect estimation from coloc-validated SNPs
- causal-genomics/fine-mapping - SuSiE / FINEMAP / CAVIAR credible sets feeding coloc.susie; LD construction protocol cross-ref
- causal-genomics/mediation-analysis - Downstream causal mediation given coloc shared causal variants
- causal-genomics/pleiotropy-detection - Distinguishing horizontal pleiotropy from shared causality
- causal-genomics/transcriptome-wide-association - TWAS / PrediXcan / FOCUS gene-level prioritization complementary to coloc
- causal-genomics/proteome-mr-drug-target - pQTL coloc + MR for drug-target prioritization
- causal-genomics/effector-gene-prioritization - Locus-to-gene with coloc, ABC, V2G integration
- population-genetics/association-testing - GWAS summary statistic generation and locus extraction
- population-genetics/linkage-disequilibrium - LD reference panels for coloc.susie and PWCoCo
- variant-calling/variant-annotation - Functional annotation for variant-specific priors
- variant-calling/filtering-best-practices - Pre-coloc QC for summary stats
- differential-expression/deseq2-basics - Generating eQTL / molecular QTL counts
- single-cell/scatac-analysis - Per-cell-type chromatin context for coloc interpretation
- workflows/gwas-pipeline - Upstream GWAS analysis producing coloc input
- data-visualization/ggplot2-fundamentals - Regional and LocusCompare plot construction
