---
name: bio-causal-genomics-proteome-mr-drug-target
description: Runs cis-pQTL Mendelian randomization for drug-target validation using UKB-PPP (Olink), deCODE (SomaScan), Fenland, INTERVAL, ARIC, and FinnGen-PPP proteomes plus colocalization triangulation, phenome-wide on-target adverse-effect scans, cross-platform Olink/SomaScan replication, and PAV (protein-altering variant) sensitivity. Use when nominating or de-risking a drug target from plasma-proteome GWAS, mimicking pharmacological inhibition via cis-pQTL instruments, separating shared-causal from LD-confounded signal under the Schmidt 2020 cis-MR framework, screening on-target adverse phenotypes pheWAS-style, or producing publication-grade STROBE-MR plus PP.H4 evidence for a target gene.
tool_type: mixed
primary_tool: TwoSampleMR
---

## Version Compatibility

Reference examples tested with: TwoSampleMR 0.5.11+, MendelianRandomization 0.10+, MR-PRESSO 1.0+, coloc 5.2.3+, susieR 0.12.35+, ieugwasr 1.0+, plink2 2.00a5+, R 4.4+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `plink2 --version`; VEP `vep --help`

If code throws OAuth or rate-limit errors from OpenGWAS, or a missing `dataset$N` from coloc, introspect the installed API and adapt the example rather than retrying. UKB-PPP, deCODE, and Fenland summary statistics changed file layouts between 2023 and 2025; verify column headers before passing into `format_data()`.

# Proteome-Wide Drug-Target Mendelian Randomization

**"Does genetically lowering plasma protein X cause a change in disease Y, mimicking a drug?"** -> Use cis-pQTLs in the gene window for protein X as instruments under the Schmidt 2020 framework (Nat Commun 11:3255), restrict the exclusion-restriction violation to the geometric neighbourhood of the encoding gene, triangulate with colocalization (3-tier PP.H4 ladder below) and cross-platform replication (Olink vs SomaScan), and flag protein-altering-variant (PAV) confounding. A single significant cis-MR estimate is necessary but not sufficient for a drug-target claim; the operational bar is MR + coloc + cross-platform agreement + PAV-excluded sensitivity.

### PP.H4 Three-Tier Threshold Ladder

| Tier | PP.H4 | Use case |
|------|-------|----------|
| Suggestive | >= 0.7 | Open Targets / exploratory; consistent with shared-causal |
| Standard publication | >= 0.8 | Wallace 2020 PLoS Genet 16:e1008720; most peer-reviewed pubs |
| Industry / clinical | >= 0.95 | Drug-claim grade; pharma internal target-validation standard |

Operational rule: drug-target nomination requires PP.H4 >= 0.8 minimum; industry-grade clinical claim requires PP.H4 >= 0.95 plus the full triangulation panel.

- R (canonical): `TwoSampleMR::mr()` orchestrates the cis-IVW + Egger + median + Wald-ratio panel
- R (correlated cis-pQTLs in a window): `MendelianRandomization::mr_ivw(model='default', correl=TRUE, correl.x=ld_matrix)`
- R (triangulation): `coloc::coloc.abf()` or `coloc::coloc.susie()` on the same cis-window
- pheWAS: `ieugwasr::associations()` against the OpenGWAS catalogue, looped over outcomes
- VEP CLI: annotate every cis-pQTL with `vep --species homo_sapiens --canonical --check_existing` for PAV flagging

## Data Source Taxonomy

| pQTL dataset | Platform | N | Proteins | Reference | Fails when |
|--------------|----------|---|----------|-----------|------------|
| UKB-PPP | Olink Explore 3072 (antibody PEA) | 54,219 | 2923 (54k primary cis-pQTLs across studies; 14,287 primary pQTLs across cis+trans) | Sun 2023 Nature 622:329 | Target not on Olink panel; ancestry mostly EUR; antibody epitope may miss isoforms |
| deCODE | SomaScan v4 (aptamer SOMAmer) | 35,559 | 4907 | Ferkingstad 2021 Nat Genet 53:1712 | Population isolate; LD differs from outbred EUR; aptamers can be PAV-confounded |
| Fenland | SomaScan v4 | 10,708 | 4775 | Pietzner 2021 Science 374:eabj1541 | UK Fenland-specific ascertainment; SomaScan caveats |
| INTERVAL | SomaScan v3 (older) | 3301 | 3622 | Sun 2018 Nature 558:73 | Smaller N; older SomaScan version; useful only as replication |
| ARIC | SomaScan v4 | ~7000 (multi-ancestry sub-cohorts) | 4877 | Zhang 2022 Nat Genet 54:593 | Stratify by ancestry; do not pool |
| FinnGen-PPP | Olink Explore | ~12,000 | ~3000 | FinnGen DF12 (2024 release) | Finnish-specific allele frequencies; not always meta-analyse with UKB |
| AGES-Reykjavik | SomaScan v4 | ~5400 | 4782 | Emilsson 2018 Science 361:769 | Elderly Icelandic cohort; ascertainment bias |
| Olink Explore 1536 disease-specific cohorts (CARDIoGRAMplusC4D, etc) | Olink Explore subset | varies | varies | per-study | Lower N per cohort; use as replication |

Methodology evolves; check the UKB-PPP portal (ukb-ppp.gwas.eu), deCODE Genetics summary-stat releases, and the eQTL Catalogue / GTEx Portal for the current data version. UKB-PPP was substantially re-released in 2024 (extended ancestry meta-analyses); pin a download date in the methods.

### Platform and Cohort Versioning (2024-2026)

- **UKB-PPP Phase 2 (2024)** -- cross-ancestry meta-analyses; ~54k EUR plus multi-ancestry expansion; file layouts differ from Phase 1 (per-protein parquet vs flat TSV); the 14,287 primary pQTL count from Phase 1 shifts in Phase 2 results. Verify the freeze used.
- **FinnGen-PPP** -- Olink Explore platform; DF12 (2024) release is current; Finnish-specific allele frequencies require ancestry-aware downstream analysis.
- **AoU pQTL (2024)** -- All-of-Us proteomics; pre-symptomatic-cohort design strength for reverse-causation control and longitudinal follow-up.
- **SomaScan v4 vs v5** -- v5 expands to ~11k SOMAmers (vs ~5k in v4); binding consistency for shared SOMAmers documented in deCODE / SomaLogic technical notes but not guaranteed; pin platform version.
- **Olink Explore HT** -- ~5,400 proteins (vs ~3,072 in Explore Expansion, ~1,536 in Explore 1536); UKB-PPP Phase 1 used Explore 3072 -- do not assume Explore HT coverage when reading Phase 1 papers.
- Pin specific platform version + release date in the methods section of every cis-MR drug-target manuscript.

## Cis-MR Methodological Taxonomy

| Method | Cis-window assumption | Min cis-pQTLs | Strength | Fails when |
|--------|------------------------|----------------|----------|------------|
| Single cis-pQTL Wald ratio | Single sentinel SNP within +/-500 kb | 1 | Simplest, transparent point estimate `beta_Y/beta_X` and ratio SE | Confounded by LD-linked eQTL/pQTL of neighbour gene; no heterogeneity test |
| Cis-IVW (clumped r2 < 0.1) | Multiple weakly-correlated cis-pQTLs | 2 | Pools information, increases precision (Schmidt 2020) | r2 between pQTLs > 0.1 inflates SE under independence assumption |
| Cis-IVW with correlation `correl=TRUE` | Cis-pQTLs in moderate LD; supply LD matrix | 2 | Correct SE under correlated instruments (Burgess, Zuber, Valdes-Marquez, Sun, Hopewell 2017 Genet Epidemiol 41:714-725) | LD matrix mismatched to summary-stat ancestry |
| Cis-MR-Egger | Directional pleiotropy across cis-pQTLs | 3+ (>=10 for power) | Sensitivity for in-window directional pleiotropy | Underpowered <10 cis-pQTLs; NOME violation `I^2_GX < 0.9` |
| Cis-weighted-median | Up to 50% invalid cis-pQTLs | 3+ | Robust to a minority of bad instruments | >50% invalid cis-pQTLs |
| MR-PRESSO in cis-window | Outlier cis-pQTLs from LD-confounded neighbours | 4+ | Removes neighbour-eQTL-tagged cis-pQTLs; distortion test (Verbanck 2018) | <4 instruments; underpowered global test |
| Robust cis-MR (Patel 2023 Biometrics 79:3458) | Accommodates LD between cis-pQTLs without full LD matrix | 2+ | Methodologically modern alternative when LD reference is uncertain | Newer; benchmarks evolving |
| Generalized cis-IVW with correlated SNPs | Joint multivariable cis-window | 2+ | Sound under any LD provided matrix supplied | Numerical instability when r2 ~ 1 (collinear) |
| coloc.susie + Wald-ratio per CS | Multiple independent cis-signals (allelic heterogeneity) | 1+ per credible set | Per-signal MR + per-signal PP.H4 | Requires ancestry-matched LD; spurious CS under mismatch |

Verify against Burgess 2023 *Wellcome Open Res* "Guidelines for performing Mendelian randomization" (v3+) and the Open Targets Genetics drug-target pipeline (Mountjoy 2021) before pinning a primary method for a publication.

## Decision Tree by Scenario

| Scenario | Primary method | Triangulation | Why |
|----------|----------------|----------------|-----|
| Single drug target, single sentinel cis-pQTL, single outcome | Wald ratio | coloc.abf PP.H4 + cross-platform replication + PAV flag | Minimum publishable cis-MR; simplest and most transparent |
| Single drug target, multiple independent cis-pQTLs, single outcome | Cis-IVW correl=TRUE with in-window LD | coloc.susie per credible set + cross-platform | Pools signal; correctly handles within-window LD |
| Drug target with secondary independent cis-signal (allelic heterogeneity) | coloc.susie + per-CS Wald ratio | Compare effect direction across CSs | Each independent signal is its own instrument; report each |
| Phenome-wide MR on a single target | Loop Wald ratio or cis-IVW across hundreds of outcomes | Bonferroni over outcomes; coloc PP.H4 on top hits | On-target adverse-effect discovery (e.g. PCSK9 -> T2D) |
| On-target adverse-effect scan (already-marketed drug) | pheWAS cis-MR vs all FinnGen / OpenGWAS phenotypes | Bonferroni + coloc + clinical-event registry | Re-derives known and novel on-target effects |
| Cross-platform replication required for clinical claim | Run independently on UKB-PPP (Olink) AND deCODE (SomaScan) | Direction agreement + magnitude within 2x | Single platform never sufficient for therapeutic decision |
| Trans-pQTL "wants" to be an instrument | Refuse | -- | Trans = horizontal pleiotropy by definition; use only as confirmatory |
| Target gene not on Olink panel | deCODE / Fenland SomaScan only | Cross-replicate across two SomaScan cohorts | Olink panel is gated; do not infer "untested" as null |
| Target in cis-region with strong neighbour eQTL | coloc.susie + coloc to non-target gene's pQTL/eQTL | Drop cis-pQTLs that coloc with non-target | Mandatory: must rule out neighbour-gene mediation |
| Sample overlap (UKB-PPP exposure + UKB phenotype outcome) | MR-RAPS with one-sample-equivalent correction OR independent outcome (FinnGen, BBJ) | Repeat in non-overlapping cohort | Pretending overlap is two-sample inflates effect estimates |

## Per-Method Failure Modes

### Olink vs SomaScan platform discordance

**Trigger:** A cis-pQTL effect size or even direction differs between UKB-PPP (Olink antibody) and deCODE/Fenland (SomaScan aptamer) for the same protein.

**Mechanism:** Olink uses paired antibody proximity-extension assay (PEA) binding distinct epitopes; SomaScan uses single-aptamer SOMAmer binding a folded epitope. A missense SNP that alters one epitope produces an apparent pQTL on that platform only (a pseudo-PAVQTL / aptamer-affinity QTL, AAVQTL). Splice/isoform differences also produce platform-specific signal. Pietzner 2021 and Sun 2023 reported ~15-30% protein-level discordance.

**Symptom:** Cis-MR significant on one platform, null on the other; or significant on both but opposite direction.

**Fix:** Replicate every cis-MR claim on at least one alternate platform; annotate cis-pQTLs with VEP and flag missense / nonsense / splice variants in the gene; consult Olink and SomaScan documentation for the protein's antibody / aptamer binding region; perform a PAV-excluded sensitivity analysis and report both estimates. If platforms disagree irreconcilably, report the protein as platform-discordant and do NOT advance for clinical claim.

**Olink panel coverage pre-check:** Olink Explore 3072 covers ~3,000 proteins; Explore Expansion ~3,000; Explore HT ~5,400. Always confirm the target is on the panel BEFORE running cis-MR: `grep -i <GENE> olink_panel_proteins.tsv` (panel manifest from olink.com). Sun 2023 Nature supplementary Table S1 lists every UKB-PPP-covered protein explicitly; absence from Table S1 means the target was not measured in Phase 1 (regardless of biology) and the analysis must use SomaScan.

### LD-based pleiotropy in the cis-window

**Trigger:** A cis-pQTL for target gene X is also a strong eQTL or pQTL for a neighbouring gene within +/-500 kb.

**Mechanism:** The instrument's effect on outcome may be mediated by the neighbour gene's protein, not by target X. Cis-window proximity does NOT guarantee specificity.

**Symptom:** Colocalization with the non-target gene's pQTL or eQTL returns PP.H4 >= 0.7; cis-MR using only "clean" cis-pQTLs (those not coloc'd with neighbours) gives a substantially different estimate.

**Fix:** For every cis-pQTL, run colocalization against all eQTL/pQTL signals within +/-500 kb in the relevant tissue; drop cis-pQTLs that coloc (PP.H4 >= 0.5) with any non-target gene; coloc.susie is preferred when multiple credible sets exist in the window. Reports must list which cis-pQTLs were retained and the rationale.

### Reverse causation from disease state on plasma protein

**Trigger:** The outcome trait elevates the protein as a downstream consequence (e.g. CRP elevated in coronary disease patients; TNF in autoimmune disease).

**Mechanism:** Plasma proteomes measured in observational cohorts (including UKB-PPP) include people who have or will develop the outcome; the observed protein-disease association may be downstream not upstream.

**Symptom:** Observational protein-disease association is large; cis-MR estimate is much smaller, null, or in the opposite direction.

**Fix:** Apply Steiger filtering on each cis-pQTL (`steiger_filtering()`); restrict to pre-symptomatic samples where possible (pediatric or early-adult cohorts); replicate in longitudinal cohorts measuring protein years before disease onset. Note: Steiger has its own caveat under unmeasured confounding (Hemani Tilling 2022 Wellcome Open Res 7:14); cross-validate via bidirectional cis-MR.

```r
library(TwoSampleMR)
dat <- harmonise_data(exposure_pQTL, outcome_GWAS)
dat <- steiger_filtering(dat)
dat_forward <- dat[dat$steiger_dir, ]   # drop reverse-direction SNPs
dir_test <- directionality_test(dat_forward)
```

### PAV (protein-altering-variant) confound

**Trigger:** A cis-pQTL is itself a missense, nonsense, frameshift, splice-site, or stop-gain variant in the target gene.

**Mechanism:** The variant changes the protein sequence, which can change antibody affinity (Olink) or SOMAmer affinity (SomaScan) without changing the actual protein abundance in plasma. The pQTL appears strong but reflects measurement artifact.

**Symptom:** The strongest cis-pQTL is a coding variant; cis-MR effect magnitude shrinks substantially when PAVs are excluded; the pQTL is platform-specific.

**Fix:** Annotate ALL cis-pQTLs with Ensembl VEP (`--check_existing --canonical`). Tabulate every cis-pQTL's most-severe consequence. Run two cis-MR analyses: (a) all cis-pQTLs, (b) PAV-excluded. Report both; require concordance for a publishable claim (Sun 2023 supplementary). For aptamer panels, also annotate the SOMAmer binding-region overlap if available.

**PAV-excluded concordance rule:** Concordance between all-cis and PAV-excluded estimates requires (i) effect direction preserved, (ii) |effect| within 2x of the all-cis estimate, AND (iii) p-value still nominally significant (P < 0.05) after PAV exclusion. If ANY of the three criteria fails, report both estimates and downgrade the claim from "drug-target" to "suggestive cis association requiring orthogonal confirmation."

### Trans-pQTL pleiotropy if used as instrument

**Trigger:** Including trans-pQTLs (outside +/-500 kb of the gene) in the instrument set to gain power.

**Mechanism:** A trans-pQTL acts through some other gene's protein that then regulates target X; using it as an instrument violates the exclusion-restriction by definition (horizontal pleiotropy).

**Symptom:** Effect estimate shifts when trans-pQTLs are added; Egger intercept significant.

**Fix:** Restrict instrument set to cis-window only (Schmidt 2020). Trans-pQTLs may be confirmatory ("does the regulator gene also predict outcome?") but never primary instruments for a drug-target claim.

### Sample overlap when both ends are UKB

**Trigger:** Using UKB-PPP for the exposure (Olink protein) AND a UKB phenotype (HES, cancer registry, ICD-10) for the outcome.

**Mechanism:** The same individuals contribute to both summary statistics; weak-instrument bias is now one-sample-equivalent and points TOWARD the confounded observational estimate, not the null.

**Symptom:** UKB-on-UKB cis-MR effect is substantially larger than the same protein-disease pair tested in an independent cohort.

**Fix:** Use an independent outcome cohort whenever possible (FinnGen, Biobank Japan, MVP). If UKB-on-UKB is necessary, apply MR-RAPS with one-sample-equivalent treatment OR the Burgess 2016 (Genet Epidemiol 40:597) sample-overlap correction. Document the overlap fraction. MRlap (Mounier 2023 Genet Epidemiol 47:314) is the recommended modern tool for joint sample-overlap and winner's-curse correction; see causal-genomics/mendelian-randomization for the canonical implementation.

## Triangulation Requirement (Operational Postdoc Rule)

A drug-target causal claim that survives peer review and informs pharmacology requires MULTIPLE concordant streams, not a single significant cis-MR p-value:

1. **Cis-MR estimate** with `P < 0.05 / N_proteins` (Bonferroni proteome-wide ~ 1.7e-5 for 2923 Olink proteins) OR `P < 0.05 / N_outcomes` (Bonferroni pheWAS-wide) -- depending on the testing regime
2. **Colocalization PP.H4 >= 0.8** between cis-pQTL and outcome GWAS at the gene locus -- standard publication tier; >= 0.95 for industry-grade claim (cross-reference causal-genomics/colocalization-analysis)
3. **Cross-platform replication** -- significant on both Olink (UKB-PPP) AND SomaScan (deCODE or Fenland); direction agreement is mandatory, magnitude within 2x is acceptable
4. **Cross-cohort replication** (ideal but not strictly required) -- UKB-PPP -> deCODE -> FinnGen-PPP step-up
5. **PAV-excluded sensitivity** -- the cis-MR survives when missense / nonsense / splice / aptamer-binding-region SNPs are dropped
6. **No neighbour-gene coloc** -- no cis-pQTL in the instrument set coloc-shares a causal variant with any non-target gene's eQTL or pQTL within +/-500 kb
7. **Open Targets L2G concordance** -- Open Targets Platform locus-to-gene score for the target gene >= 0.5 at the disease GWAS lead (cross-reference causal-genomics/effector-gene-prioritization)

Operational claim ladder: cis-MR significant alone = exploratory; +coloc PP.H4 >= 0.7 = consistent with shared-causal; +cross-platform = consistent across detection chemistries; +PAV-excluded + neighbour-clear = publication-grade target nomination; +cohort replication + Open Targets L2G >= 0.5 = clinical-pharmacology-grade. Clinical-pharmacology claims require ALL 6 original criteria PLUS L2G concordance.

## Phenome-Wide Drug-Target MR

**Goal:** For a single drug target (single protein), test causal effect across hundreds of outcomes to discover on-target adverse effects.

**Approach:** Hold cis-pQTL instrument set fixed; loop outcome over OpenGWAS catalogue or FinnGen DF12; multi-test correct over outcomes.

```r
library(TwoSampleMR); library(ieugwasr)

target_pqtl <- read.table('pcsk9_cis_pqtls.tsv', header = TRUE)
exposure_dat <- format_data(target_pqtl, type = 'exposure',
    snp_col = 'SNP', beta_col = 'BETA', se_col = 'SE',
    effect_allele_col = 'A1', other_allele_col = 'A2', eaf_col = 'EAF', pval_col = 'P')

curated_endpoints <- read.table('finngen_DF12_endpoints.tsv', header = TRUE)  # ~3000 curated endpoints from finngen.fi
outcomes <- available_outcomes()
outcomes_filt <- subset(outcomes, id %in% curated_endpoints$id & sample_size >= 50000 & population == 'European')

results <- lapply(outcomes_filt$id, function(out_id) {
    outcome_dat <- extract_outcome_data(snps = exposure_dat$SNP, outcomes = out_id)
    if (nrow(outcome_dat) < 2) return(NULL)
    dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
    mr(dat, method_list = 'mr_ivw')
})

results_df <- do.call(rbind, Filter(Negate(is.null), results))
n_tests <- nrow(curated_endpoints)
results_df$p_bonf <- pmin(results_df$pval * n_tests, 1)
top_hits <- subset(results_df, pval < 0.05 / n_tests)   # 0.05 / 3000 = 1.7e-5
```

Document the curated endpoint list (FinnGen DF12, Open Targets curated trait map, or a manuscript-specific phecode hierarchy) in methods. The `outcomes_filt$id[1:200]` pattern is a debug shortcut, not a defensible pheWAS protocol. The PCSK9 -> T2D signal (Schmidt 2017 Lancet Diabetes Endocrinol 5:97) was discovered exactly via curated-endpoint pheWAS: cis-MR of LDL-lowering instruments revealed on-target T2D risk before clinical trials confirmed it. Drug-target pheWAS is the canonical use case.

## Cis-MR Standard Workflow

**Goal:** Produce a defensible cis-MR estimate with triangulation, given a target gene and an outcome.

**Approach:** Extract cis-pQTLs in +/-500 kb of the gene -> compute F per instrument from exposure -> clump within window at r2 < 0.1 -> harmonise with outcome -> run TwoSampleMR -> run coloc.abf on the same window -> PAV-annotate via VEP -> report panel.

```r
library(TwoSampleMR); library(coloc); library(ieugwasr)

cis_window_kb <- 500  # +/- 500 kb per Schmidt 2020 standard cis-window
gene_chr <- 1; gene_start <- 55039548; gene_end <- 55064852  # PCSK9 hg38

pqtl <- read.table('ukbppp_pcsk9.tsv', header = TRUE)
pqtl_cis <- subset(pqtl, CHR == gene_chr &
    POS > (gene_start - cis_window_kb * 1000) &
    POS < (gene_end + cis_window_kb * 1000) &
    P < 5e-8)

pqtl_cis$f_stat <- (pqtl_cis$BETA / pqtl_cis$SE)^2  # F from exposure (Burgess 2011)
pqtl_cis <- subset(pqtl_cis, f_stat >= 10)  # Staiger-Stock 1997 weak-IV floor

exposure_dat <- format_data(pqtl_cis, type = 'exposure',
    snp_col = 'SNP', beta_col = 'BETA', se_col = 'SE',
    effect_allele_col = 'A1', other_allele_col = 'A2', eaf_col = 'EAF', pval_col = 'P')

clumped <- ld_clump(
    dplyr::tibble(rsid = exposure_dat$SNP, pval = exposure_dat$pval.exposure),
    clump_r2 = 0.1, clump_kb = cis_window_kb,  # cis-MR clumping per Schmidt 2020
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = '1kg_EUR/EUR'
)
exposure_dat <- subset(exposure_dat, SNP %in% clumped$rsid)

outcome_dat <- read_outcome_data('cad_gwas.tsv', snps = exposure_dat$SNP, sep = '\t',
    snp_col = 'SNP', beta_col = 'BETA', se_col = 'SE',
    effect_allele_col = 'A1', other_allele_col = 'A2', eaf_col = 'EAF', pval_col = 'P')

dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)

primary <- mr(dat, method_list = c('mr_wald_ratio', 'mr_ivw',
                                    'mr_egger_regression', 'mr_weighted_median'))

# Triangulation step 1: colocalization on the same window
gwas_window <- read.table('cad_gwas_pcsk9_window.tsv', header = TRUE)
pqtl_window <- read.table('ukbppp_pcsk9_full_window.tsv', header = TRUE)

coloc_res <- coloc.abf(
    dataset1 = list(beta = pqtl_window$BETA, varbeta = pqtl_window$SE^2,
                    snp = pqtl_window$SNP, type = 'quant', N = 54219, sdY = 1),
    dataset2 = list(beta = gwas_window$BETA, varbeta = gwas_window$SE^2,
                    snp = gwas_window$SNP, type = 'cc', N = 122733, s = 0.34),
    p1 = 1e-4, p2 = 1e-4, p12 = 5e-6  # conservative p12 for drug-target
)

cat('PP.H4 =', coloc_res$summary['PP.H4.abf'], '\n')
```

## Cis-IVW with Correlated Instruments

**Goal:** When several cis-pQTLs in the window are correlated (r2 between 0.1 and 0.7), use the generalized IVW that takes the LD matrix as an explicit parameter.

**Approach:** Compute or load the LD matrix in the cis-window from an ancestry-matched plink reference; supply to `MendelianRandomization::mr_ivw(correl = TRUE, correl.x = ld_matrix)`.

```r
library(MendelianRandomization); library(ieugwasr)

ld <- ld_matrix(exposure_dat$SNP, bfile = '1kg_EUR/EUR', plink_bin = genetics.binaRies::get_plink_binary())

mr_obj <- mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                   by = dat$beta.outcome, byse = dat$se.outcome,
                   corr = ld)

result_correl <- mr_ivw(mr_obj, model = 'default', correl = TRUE)
```

Numerical caveat: when any pair of cis-pQTLs has r2 ~ 1 (e.g. perfect proxies), the LD matrix is rank-deficient and the SE explodes. Pre-prune at r2 < 0.95.

**API caveat for `MendelianRandomization::mr_ivw`:** Correlation between cis-pQTLs is supplied at MRInput construction via the `corr` argument (note: in older docs it appeared as `correlation=`; the canonical argument in current MendelianRandomization R is `corr=`): `mr_input(bx, bxse, by, byse, corr = ld_matrix)`. `mr_ivw()` then reads the correlation slot directly; passing `correl = TRUE` as an argument is redundant when MRInput already has a non-NA correlation matrix. A common bug is supplying both, which produces inconsistent behaviour across versions; prefer the MRInput-slot approach and treat `correl = TRUE` as a legacy flag.

### Patel 2023 Robust cis-MR with Correlated Instruments

**Goal:** Accommodate LD between cis-pQTLs when the in-sample LD reference is uncertain or when an explicit full LD matrix is not reliably available.

**Approach:** Patel, Gill, Newcombe, Burgess 2023 *Biometrics* 79:3458-3471 propose a robust cis-MR estimator that uses correlated-IV principles without requiring the full LD matrix as a known input; the estimator down-weights instruments whose correlation structure is poorly captured.

```r
library(MendelianRandomization)

mr_obj <- mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                   by = dat$beta.outcome, byse = dat$se.outcome,
                   corr = ld)
robust_res <- mr_ivw(mr_obj, model = 'default', robust = TRUE, penalized = TRUE)
```

Decision rule: standard cis-IVW when post-clumped LD r2 < 0.1; correlated-IV cis-IVW (Burgess 2017 Genet Epidemiol) when r2 between 0.1 and 0.7 AND ancestry-matched LD matrix is trustworthy; Patel 2023 robust cis-MR when LD persists after clumping AND the LD matrix may be misspecified (e.g. ancestry-mismatched reference, finite-sample noise in panel of < 500 individuals). Benchmarks for Patel 2023 vs Burgess 2017 are still evolving; report both as sensitivity when feasible.

## PAV Annotation via VEP

**Goal:** Tag every cis-pQTL with its most severe coding consequence to enable a PAV-excluded sensitivity panel.

**Approach:** Format SNPs as VEP input, run VEP, parse the `Consequence` column.

```bash
echo -e "chr1\t55039548\t.\tG\tT" > cis_pqtls.vcf
vep --species homo_sapiens --assembly GRCh38 --canonical --check_existing \
    --input_file cis_pqtls.vcf --output_file pqtl_vep.tsv --tab --force_overwrite
```

PAV consequences to flag (drop in sensitivity analysis): `missense_variant`, `stop_gained`, `stop_lost`, `frameshift_variant`, `splice_acceptor_variant`, `splice_donor_variant`, `start_lost`, `protein_altering_variant`. For aptamer panels, additionally consider `synonymous_variant` within the SOMAmer-binding region (rare but documented).

## Cis-Window Width: 500 kb vs 1 Mb Decision

- **Default ±500 kb** (Schmidt 2020 Nat Commun 11:3255 standard) -- balances cis-specificity against power; matches Open Targets Genetics defaults
- **Widen to ±1 Mb** when: (a) target gene has a documented distal regulatory element in ENCODE-rE2G or ABC enhancer-gene maps; (b) < 2 genome-wide-significant pQTLs are present in ±500 kb; (c) the target gene is unusually large (gene body > 500 kb itself, e.g. DMD, RBFOX1)
- **Narrow to ±250 kb** when high-density cis-eQTL background causes multi-gene pleiotropy concerns (e.g. HLA region; gene-dense pericentromeric loci)

Pre-specify the window in the methods section; window-width sensitivity is a recognized peer-review pushback (see Anticipated Reviewer Pushback below).

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Bonferroni for cis-MR pheWAS | Standard | `P < 0.05 / N_outcomes` for on-target adverse-effect scan |
| Bonferroni for proteome-wide cis-MR | Standard | `P < 0.05 / N_proteins` ~ 1.7e-5 for 2923 Olink proteins; ~1e-5 for 4907 SomaScan |
| Coloc PP.H4 >= 0.7 (suggestive) | Open Targets Genetics; Mountjoy 2021 Nat Genet 53:1527 | Exploratory; consistent with shared-causal |
| Coloc PP.H4 >= 0.8 (publication) | Wallace 2020 PLoS Genet 16:e1008720 | Standard peer-reviewed publication bar |
| Coloc PP.H4 >= 0.95 (industry) | Pharma internal target-validation standard | Drug-claim grade; clinical-pharmacology bar |
| Cis-pQTL F >= 10 | Staiger & Stock 1997; Burgess 2011 | Weak-instrument floor |
| Cis-window +/-500 kb | Schmidt 2020 Nat Commun 11:3255 | Standard cis definition; some pipelines use 1 Mb |
| r2 < 0.1 clumping in cis-window | Schmidt 2020 | Reduces LD-based pleiotropy while retaining power |
| N >= 2 pQTL datasets in agreement | Best-practice (UKB-PPP + deCODE) | Cross-platform replication mandatory for clinical claim |
| PAV-excluded sensitivity | Sun 2023 supplementary | Required for clinical claim; Olink/SomaScan vulnerable to PAV artifact |
| Neighbour-gene coloc PP.H4 < 0.5 | Operational | Cis-pQTL must NOT coloc with non-target gene; drop if it does |
| Steiger p > 0.05 (correct direction) | Hemani 2017 PLoS Genet 13:e1007081 | Directionality check; subject to Hemani Tilling 2022 caveat |
| Sample-overlap correction if both ends UKB | Burgess 2016 Genet Epidemiol 40:597 | One-sample-equivalent bias correction |

## Reconciliation: When Evidence Streams Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Cis-MR significant on Olink, null on SomaScan | Platform-specific epitope/aptamer artifact | PAV-annotate; flag protein as platform-discordant; do not claim |
| Cis-MR sig, coloc PP.H4 < 0.5 | LD-confounded signal; not shared causal | Downgrade; cis-MR alone insufficient; do not claim drug target |
| Coloc PP.H4 high, cis-MR null | Underpowered cis-MR (few cis-pQTLs, weak F) OR shared-causal for non-causal protein | Inspect cis-pQTL strength; consider increasing window to 1 Mb |
| Cis-MR sig with all pQTLs, null after PAV exclusion | PAV artifact dominating instrument | Report PAV-excluded as primary; original as supplementary |
| Cis-MR sig in UKB-PPP -> UKB outcome, null in FinnGen outcome | Sample-overlap one-sample bias | Treat FinnGen as truth; report UKB-on-UKB as biased upward |
| Cis-pQTL coloc with non-target gene's eQTL (PP.H4 >= 0.7) | Neighbour-gene mediation | Drop this cis-pQTL; re-run cis-MR with clean instruments |
| Two independent cis-pQTLs give opposite direction effects | Allelic heterogeneity with distinct biology | Run coloc.susie per credible set; report each signal separately |
| Wald ratio at sentinel SNP differs from cis-IVW | One outlier cis-pQTL dominates IVW | Run MR-PRESSO; check Egger intercept |

**Operational rule for publication:** Cis-IVW (or Wald ratio if N=1 SNP) + coloc.abf PP.H4 >= 0.8 + cross-platform replication (Olink and SomaScan agree in direction) + PAV-excluded sensitivity concordant + L2G >= 0.5 at disease GWAS lead = drug-target nomination ready. Industry-grade clinical claim requires PP.H4 >= 0.95. Any single missing leg downgrades the claim to "consistent with" rather than "evidence for."

## Drug Repurposing and Target Nomination

The Open Targets Drug platform (Ochoa 2021 Nucleic Acids Res 49:D1302) integrates approved-drug-target relationships, cis-pQTL/cis-eQTL MR, and locus-to-gene (L2G) scores (Mountjoy 2021). A target is nominated for repurposing when:

- L2G score at the GWAS lead points to the target gene
- Cis-MR estimate concordant with disease direction (lower protein -> lower disease for inhibitor candidate)
- Coloc PP.H4 >= 0.7 between target's cis-pQTL and disease GWAS
- A licensed drug exists that modulates the target
- The on-target adverse-effect pheWAS is acceptable

The PCSK9 monoclonal-antibody story is the canonical positive example; the CETP-inhibitor story (Joshi 2020) is a canonical cautionary tale (cis-MR underestimated trial result due to off-target effects).

Cross-validate target nominations against the Comparative Toxicogenomics Database (CTD; ctdbase.org) and the Drug-Gene Interaction Database (DGIdb; dgidb.org) for established drug-target evidence and tractability annotations. STROBE-MR (Skrivankova 2021 JAMA 326:1614) provides the 20-item checklist for drug-target MR reporting; see causal-genomics/pleiotropy-detection for the full table.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `harmonise_data` drops most SNPs | EAF columns missing or palindromic at MAF~0.5 | Provide EAF; use `action = 2` (default) or `action = 3` for strictest |
| F-statistic from outcome | Computed `beta.outcome / se.outcome` | F must come from exposure (Burgess 2011) |
| Trans-pQTL included as instrument | Filtered only by P, not by genomic position | Restrict to +/-500 kb of target gene |
| PAV cis-pQTL artifact | Did not VEP-annotate | Run VEP; flag missense/nonsense/splice; report PAV-excluded sensitivity |
| Neighbour-gene mediation | Did not coloc cis-pQTL against neighbours | Run coloc.susie or coloc.abf vs non-target eQTL/pQTL in window |
| OpenGWAS OAuth failure | Token expired (OpenGWAS auth tightened 2024) | Use local plink + 1KG bfile via `ieugwasr::ld_clump(..., bfile=...)` |
| `mr_ivw(correl=TRUE)` SE explodes | Two cis-pQTLs in near-perfect LD | Pre-prune at r2 < 0.95; or drop redundant proxy |
| Cis-IVW null but Wald ratio at lead significant | Inclusion of weak / outlier cis-pQTLs | Tighten clumping; run MR-PRESSO outlier test |
| Sample overlap unreported | Both GWAS from UKB; analyst assumed two-sample | Apply Burgess 2016 correction OR use MR-RAPS OR move outcome to FinnGen |
| Coloc PP.H3 dominant | Multiple causal in moderate LD | Switch to coloc.susie with ancestry-matched LD |
| Olink and SomaScan disagree | Platform-specific epitope artifact | Annotate PAV; flag as platform-discordant; do not claim drug-target |

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Cross-platform replication?" | Olink (UKB-PPP) and SomaScan (deCODE / Fenland) agreement reported; direction match mandatory, magnitude within 2x; PAV-excluded sensitivity included |
| "PAV check?" | Ensembl VEP annotation of every cis-pQTL; PAV-excluded sensitivity reported alongside primary; concordance criteria (direction + 2x magnitude + nominal P) documented |
| "Neighbour-gene coloc?" | coloc.abf of each cis-pQTL against all eQTLs and pQTLs within ±500 kb in the relevant tissue; PP.H4 < 0.5 for off-target genes required for instrument retention |
| "Reverse causation?" | Steiger filter applied per cis-pQTL; bidirectional cis-MR run; pre-symptomatic-cohort replication (AoU, longitudinal sub-studies) reported when available |
| "Sample overlap (UKB-on-UKB)?" | MRlap correction applied (Mounier 2023); or outcome moved to independent cohort (FinnGen DF12, MVP, Biobank Japan); overlap fraction documented |
| "Industry-grade threshold?" | PP.H4 >= 0.95 reported for drug-claim grade; >= 0.8 for standard publication; 3-tier ladder pre-specified in methods |
| "OT L2G concordance?" | Open Targets Platform L2G score >= 0.5 at the disease GWAS lead reported (cross-reference effector-gene-prioritization) |
| "Cis-window width?" | ±500 kb default (Schmidt 2020) pre-specified; widened to ±1 Mb only with documented distal regulatory rationale |
| "Patel 2023 vs Burgess 2017?" | Both reported when post-clumped LD r2 > 0.1; Patel 2023 robust estimator preferred when LD reference is ancestry-mismatched |

## Tool Installation Notes

```r
install.packages(c('remotes', 'MendelianRandomization', 'coloc', 'susieR', 'dplyr'))
remotes::install_github('MRCIEU/TwoSampleMR')
remotes::install_github('MRCIEU/ieugwasr')
remotes::install_github('MRCIEU/genetics.binaRies')   # bundles plink binary
remotes::install_github('rondolab/MR-PRESSO')
```

```bash
# Ensembl VEP for PAV annotation
conda install -c bioconda ensembl-vep
vep_install -a cf -s homo_sapiens -y GRCh38 -c $HOME/.vep

# 1000 Genomes EUR plink reference for local clumping / LD matrix
# Prebuilt at https://mrcieu.github.io/ieugwasr/
```

pQTL data acquisition: UKB-PPP via the UK Biobank pre-published portal (ukb-ppp.gwas.eu); deCODE via the deCODE Genetics summary-stat website with a data-use agreement; Fenland via the EBI GWAS catalog and Pietzner 2021 supplementary; OpenGWAS hosts many pre-formatted pQTL studies but always verify the upstream reference and download date.

## References

- Schmidt AF et al 2020 Nat Commun 11:3255 (cis-MR drug-target framework)
- Sun BB et al 2023 Nature 622:329 (UKB-PPP Olink Explore)
- Ferkingstad E et al 2021 Nat Genet 53:1712 (deCODE SomaScan)
- Pietzner M et al 2021 Science 374:eabj1541 (Fenland SomaScan)
- Sun BB et al 2018 Nature 558:73 (INTERVAL SomaScan)
- Zhang J et al 2022 Nat Genet 54:593 (ARIC multi-ancestry pQTL)
- Emilsson V et al 2018 Science 361:769 (AGES-Reykjavik)
- Schmidt AF et al 2017 Lancet Diabetes Endocrinol 5:97 (PCSK9 -> T2D pheWAS exemplar)
- Burgess S, Zuber V, Valdes-Marquez E, Sun BB, Hopewell JC 2017 Genet Epidemiol 41:714-725 (correlated-IV IVW)
- Burgess S et al 2016 Genet Epidemiol 40:597 (sample-overlap correction)
- Mountjoy E et al 2021 Nat Genet 53:1527 (Open Targets Genetics L2G + coloc)
- Ochoa D et al 2021 Nucleic Acids Res 49:D1302 (Open Targets Drug platform)
- Hemani G & Tilling K 2022 Wellcome Open Res 7:14 (Steiger caveat)
- Joshi PK et al 2020 Eur Heart J 41:e10 (CETP cis-MR cautionary tale)
- Verbanck M et al 2018 Nat Genet 50:693 (MR-PRESSO)
- Wallace C 2020 PLoS Genet 16:e1008720 (coloc p12 sensitivity)
- Skrivankova VW et al 2021 JAMA 326:1614 (STROBE-MR)
- Patel A, Gill D, Newcombe P, Burgess S 2023 Biometrics 79:3458-3471 (robust cis-MR with correlated instruments)
- Mounier N & Kutalik Z 2023 Genet Epidemiol 47:314 (MRlap sample-overlap + winner's-curse correction)

## Related Skills

- causal-genomics/mendelian-randomization - Parent polygenic-MR framework; cis-MR is the drug-target specialization; MRlap sample-overlap correction
- causal-genomics/colocalization-analysis - Required PP.H4 triangulation for any cis-MR drug-target claim
- causal-genomics/fine-mapping - Credible-set construction prior to coloc.susie at the cis-locus
- causal-genomics/pleiotropy-detection - MR-PRESSO / Egger diagnostics adapted to cis-window; STROBE-MR 20-item checklist
- causal-genomics/transcriptome-wide-association - eQTL-based parallel evidence for the same target
- causal-genomics/mediation-analysis - Step from cis-MR to downstream mediator pathway
- causal-genomics/effector-gene-prioritization - Open Targets L2G concordance leg of the triangulation panel
- population-genetics/association-testing - Source GWAS pipelines for pQTL discovery
- population-genetics/linkage-disequilibrium - LD-matrix construction for cis-IVW correl=TRUE
- variant-calling/variant-annotation - VEP PAV annotation for sensitivity analysis
- clinical-databases/clinvar-lookup - Pathogenic-variant context for nominated targets
