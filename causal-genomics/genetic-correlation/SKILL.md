---
name: bio-causal-genomics-genetic-correlation
description: Estimate bivariate genetic correlation (rg) between traits from GWAS summary statistics or individual-level genotypes using cross-trait LDSC, HDL, LAVA, rho-HESS, GREML-bivariate, Popcorn, and HDL-L. Use when quantifying shared genetic architecture between two traits, screening MR validity before causal inference, distinguishing global from locus-level rg, estimating trans-ancestry rg, separating partial from full causation via LCV gcp, or producing a STROBE-MR-compliant cross-trait sensitivity battery. Cross-trait LDSC intercept absorbs sample overlap and is NOT a bias; HDL is biased under sample overlap above ~5%. High rg between exposure and outcome motivates CHP-aware MR sensitivity (CAUSE, LHC-MR).
tool_type: mixed
primary_tool: ldsc
---

## Version Compatibility

Reference examples tested with: LDSC v1.0.1+ (Python 3; prefer `abdenlab/ldsc-python3` v2.0.0 -- `belowlab/ldsc` v3.0.1 README states the `--h2 / --rg / --h2-cts` CLI is broken; use Docker `jtb114/ldsc:latest` for the belowlab fallback; original `bulik/ldsc` is Python 2.7 unmaintained since 2019), HDL 1.4.0+ (R; GitHub `zhenin/HDL`), LAVA 0.1.0+ (R; GitHub `josefin-werme/LAVA`), HESS 0.5.4+ (Python; huwenboshi/hess), Popcorn 1.0+ / Popcorn-2 (Python; brielin/Popcorn), GCTA 1.94+ (GREML-bivariate), baselineLD_v2.2 / eur_w_ld_chr LD-score panels from alkesgroup.broadinstitute.org/LDSCORE, UKB-array SVD eigen reference for HDL.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `python -c 'import <module>; help(<module>)'`
- R: `packageVersion('<pkg>')` then `?function_name`
- CLI: `<tool> --version` then `<tool> --help`

If code throws an LD-score "category not found" error, an HDL reference-panel mismatch, or a LAVA locus-ID lookup failure, introspect the installed LD-score column headers and the supplied partitioning file rather than retrying with default flags.

# Genetic Correlation

**"Estimate the genetic correlation between two traits from GWAS summary statistics"** -> Decompose the bivariate genetic architecture into a single global rg (cross-trait LDSC, HDL), per-locus local rg (LAVA, rho-HESS, HDL-L), or cross-population rg (Popcorn / Popcorn-2). Genetic correlation is the central cross-trait statistic in causal genomics: it quantifies shared etiology, motivates CHP-aware MR sensitivity when high, gates LCV's gcp partial-causation parameter, and feeds into multi-trait analysis frameworks (MTAG, GenomicSEM). Tool choice is a decision about the **regime** (sumstats vs individual-level; global vs local; same-ancestry vs trans-ancestry) and the **sample-overlap structure** between input GWAS.

- CLI (LDSC, robust to overlap): `ldsc.py --rg trait1.sumstats.gz,trait2.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out rg`
- R (HDL, lower variance, requires independent samples): `HDL.rg(gwas1.df, gwas2.df, LD.path = 'UKB_array_SVD_eigen99_extraction', N0 = 0)`
- R (LAVA, local rg per locus): `process.input() -> run.univ() -> run.bivar(input, locus_id)` over ~2495 LDetect-derived loci
- CLI (rho-HESS, locus-level): `hess.py --local-rhog t1.sumstats.gz,t2.sumstats.gz --bfile <ref> --partition-file <part>.bed --chrom <chr>`
- CLI (Popcorn, trans-ancestry): `popcorn fit -v 1 --cfile cross_pop_scores.txt --sfile1 pop1.txt --sfile2 pop2.txt out`

## Algorithmic Taxonomy

| Method | Model | Input | Output | Strength | Fails when |
|--------|-------|-------|--------|----------|------------|
| Cross-trait LDSC (Bulik-Sullivan 2015 Nat Genet 47:1236) | Bivariate LD-score regression; off-diagonal absorbs rg, intercept absorbs sample overlap | Sumstats + ancestry-matched LD scores | rg, SE, intercept (overlap proxy) | Robust to sample overlap (intercept absorbs it without biasing rg); fast; calibrated EUR | Mean chi-square < 1.02 in either trait (underpowered); non-EUR sumstats with EUR LD scores |
| HDL (Ning 2020 Nat Genet 52:859) | High-Definition Likelihood; eigen-decomposition of full LD with closed-form variance | Sumstats + UKB-array SVD eigen reference (EUR N=336k) | rg, SE | ~60% lower variance than LDSC; equivalent to ~2.5x sample size; preferred when both GWAS truly independent | Sample overlap > 5% biases likelihood; only public reference panel is EUR UKB-array |
| LAVA (Werme 2022 Nat Genet 54:274) | Semi-parametric local genetic correlation per locus; PC-projected SNP effects under a local null | Sumstats + LDetect partitioning (~2495 loci) | Per-locus univariate h2 + bivariate rg + p-value | Detects heterogeneous rg masked by global cancellation; conditional + partial rg supported | Locus has too few SNPs (< 50) or low local h2 (univariate p > 0.05 in either trait); LD reference mismatch |
| rho-HESS (Shi 2017 AJHG 101:737) | Quadratic form on LD-projected effect estimates per locus | Sumstats + LDetect partition + LD reference | Per-locus rho_g + bivariate local rg | Earliest locus-level rg method; complements LAVA | Locus < 1000 SNPs; LD ref must match in-sample structure |
| HDL-L (Ning 2022 ext.) | HDL likelihood applied to local windows | Sumstats + windowed LD reference | Per-window local rg | Lower variance than rho-HESS at locus level | Same sample-overlap caveat as HDL; reference-panel coverage limited |
| GREML-bivariate (Lee 2012 Bioinformatics 28:2540) | Joint REML on bivariate GRM | Individual-level genotypes + both phenotypes | rg + SE | Gold standard at individual level; better precision than sumstats methods | Needs individual-level data on overlapping individuals OR carefully matched two-cohort; population stratification leaks |
| Popcorn (Brown 2016 AJHG 99:76) | Cross-population genetic effect (rho_ge) and impact (rho_gi) correlation under MAF-LD model | Sumstats per population + cross-population LD scores | Trans-ancestry rg + per-pop h2 | Quantifies shared causal architecture across ancestries | Effective N per population < 5000; cross-population LD score reference mismatched to GWAS ancestry |
| Popcorn-2 (Galinsky 2019 ext.) | Improved Popcorn with admixture-aware LD scores | Same as Popcorn + admixture LD scores | Same outputs, lower SE | More precise on admixed cohorts; better handling of LD mismatch | Same data-volume limit as Popcorn |
| GenomicSEM `LDSC()` (Grotzinger 2019 Nat Hum Behav 3:513) | LDSC wrapper feeding into SEM | Multiple sumstats | Genetic covariance matrix + multivariable SEM | Multi-trait extension of LDSC; common-factor and bifactor models | Same per-pair limits as LDSC; SEM identification problems |
| SUPERGNOVA (Zhang 2021 Nat Commun 12:7234) | Annotation-aware LD-block local rg with eigen-decomposed kernel | Sumstats + annotation-aware LD partition | Per-locus rg + p; orthogonal philosophy from LAVA | Different LD partitioning than LDetect; useful as triangulation against LAVA | Same chi-square floor as LDSC; non-EUR coverage limited |
| KGGSEE local rg (Zhang 2022 Bioinformatics) | Heritability-based per-region rg | Sumstats + region partition | Per-region rg | Java pipeline integrates h2 partition + local rg | Less adopted than LAVA; partition file curation overhead |

Methodology evolves; benchmark consensus shifts. Verify against the alkesgroup LDSC tutorial (current as of release), Werme 2022 LAVA paper + GitHub, and Speed 2020 *Nat Methods* model-comparison work before locking a primary method. When a claim depends on the model assumption (e.g. enrichment in shared loci), report at least two methods (e.g. LDSC global + LAVA local).

## Cross-Trait LDSC Intercept: Sample Overlap is Absorbed, Not a Bias

**The most common postdoc-level misreading:** Treating a non-zero cross-trait LDSC intercept as evidence of bias in the rg estimate.

The bivariate LDSC regression has the form:
`E[Z1 Z2] = sqrt(N1 N2) * rg * h2-product / M * LD_score + rho_overlap`

The intercept (`rho_overlap`) ABSORBS the contribution of sample overlap (correlated trait residuals on shared individuals). The slope (which carries rg) is unbiased even when overlap is non-zero. A non-zero intercept is the expected signature of sample overlap and is informative (it estimates phenotypic correlation among overlapping individuals), but it does NOT indicate that rg is contaminated.

**Operational rule:** Report the intercept alongside rg. When intercept is non-zero, document the overlap inferred (intercept = rho_phenotypic * sqrt(N_shared / N1 / N2) approximately) but do not discount rg. HDL, in contrast, assumes truly independent samples and does become biased above ~5% overlap; switch to LDSC under any non-trivial overlap.

## HDL Bias Under Sample Overlap

**The mirror image trap:** Running HDL on two GWAS that share controls or come from the same biobank.

HDL maximizes a likelihood that assumes independence of the two trait residuals after marginalizing genetics. With sample overlap, the residual correlation is non-zero and the likelihood is misspecified; bias is typically toward the rg estimate that corresponds to phenotypic correlation in the overlapping subset.

**Operational rule:** Use HDL only when sample overlap < 5%. When in doubt about overlap (e.g. two UKB-derived GWAS), compute LDSC intercept first; if intercept is materially non-zero, switch to LDSC for the primary rg estimate.

## Relationship to MR Causal Inference

Genetic correlation between an exposure and an outcome is a screening statistic, not a causal claim. High |rg| has three biological explanations:

| Explanation | Manifestation |
|-------------|----------------|
| Direct causation X -> Y | All causal SNPs of X feed through to Y; rg reflects mediated covariance |
| Shared heritable confounder (CHP) | A latent factor causes both; rg captures the shared variance with no direct edge |
| Reverse causation Y -> X | Symmetric structure; rg cannot resolve direction |

LCV's `gcp` parameter (O'Connor & Price 2018 Nat Genet 50:1728) attempts to distinguish partial from full causation: `gcp = 0` is pure correlation (no causation in tested direction), `gcp = 1` is full causation, `0 < gcp < 1` is partial causation. LCV uses the LDSC-style bivariate moments and is complementary to instrument-based MR.

**Operational rule for any MR analysis where |rg| > 0.3:** The IVW + Egger + MR-PRESSO triple is insufficient because all three are blind to CHP (Morrison 2020 Nat Genet 52:740). Add CAUSE (if sig SNPs >= 100) or LHC-MR. See causal-genomics/pleiotropy-detection for the full CHP-aware battery.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Standard EUR-EUR rg from sumstats | Cross-trait LDSC | Robust to sample overlap via intercept; standard ENCODE-equivalent default |
| Truly independent EUR samples, want maximum precision | HDL | ~60% lower variance than LDSC; equivalent to 2.5x sample size |
| Suspect heterogeneous rg across genome (e.g. neuropsychiatric pair with weak global rg) | LAVA local rg | Detects loci of shared etiology hidden by global cancellation |
| Locus-level rg with explicit LD partitioning | rho-HESS (or LAVA) | LAVA is newer and better-supported; rho-HESS remains the original framework |
| Cross-population (trans-ancestry) rg | Popcorn or Popcorn-2 | Within-population LDSC fails cross-pop; Popcorn models ancestry-specific causal architecture |
| Individual-level genotypes available | GREML-bivariate (GCTA) | Better precision than sumstats; gold standard at individual level |
| MR validity check before running MR | LDSC rg + LCV gcp | If |rg| > 0.3 add CHP-aware MR sensitivity (CAUSE / LHC-MR) |
| Multi-trait modeling (many traits jointly) | GenomicSEM `LDSC()` | SEM extension of LDSC; common-factor, bifactor, and network models |
| Sumstats with mean chi-square < 1.02 | Defer; meta-analyze to ~50k effective N first | LDSC variance explodes below this; nothing else fixes underpower |
| Confirmatory after a single LAVA hit | LDSC global + bidirectional MR + colocalization | Triangulate; LAVA flags shared etiology but does not establish causation |

## Per-Method Failure Modes

### Cross-trait LDSC intercept misread as bias

**Trigger:** Reader (collaborator, reviewer) sees a non-zero LDSC intercept and reports the rg as "biased by sample overlap".

**Mechanism:** The bivariate LDSC model partitions covariance between traits into a slope (rg-driven, scales with LD score) and an intercept (overlap-driven, constant in LD score). The slope is what carries rg, and it remains unbiased regardless of intercept value (Bulik-Sullivan 2015 Nat Genet 47:1236, Methods). Confusing the intercept with bias on the slope is a routine misinterpretation.

**Symptom:** Reviewer comment requesting "correction for sample overlap" when LDSC was already used; collaborator suggesting switching to HDL because intercept is non-zero.

**Fix:** Report rg with SE and the intercept as a separate statistic; cite Bulik-Sullivan 2015 Methods explicitly; do NOT switch to HDL (which is the wrong direction since HDL is the method that breaks under overlap, not LDSC).

### HDL bias with sample overlap

**Trigger:** Running HDL on two GWAS that share > 5% of individuals (e.g. two UKB-derived GWAS, two MVP-derived GWAS, GWAS reusing controls).

**Mechanism:** HDL likelihood assumes independent trait residuals after marginalizing genetics. Overlap induces non-zero residual correlation; the likelihood is misspecified and the estimate is pulled toward phenotypic correlation in the shared subset (Ning 2020 Nat Genet 52:859 Supplement).

**Symptom:** HDL rg differs substantially from cross-trait LDSC rg; HDL CI is narrower than expected from N alone; LDSC intercept (which absorbs overlap) is materially non-zero.

**Fix:** Compute LDSC intercept first as an overlap diagnostic; if non-zero, switch to LDSC as primary. HDL remains valid only when the two GWAS draw from non-overlapping cohorts (verify by cohort identifier, not just by file source).

### Non-EUR ancestry mismatch with EUR LD scores

**Trigger:** Running LDSC on a non-EUR GWAS (or admixed sample) with the default `eur_w_ld_chr/` reference.

**Mechanism:** LD scores are ancestry-specific; mean LD per SNP differs across populations and bivariate moments use the wrong null. h2 and rg estimates are systematically biased; intercept can be inflated.

**Symptom:** LDSC ratio is unusually high (>0.2); per-chromosome estimates wildly heterogeneous; total h2 mismatches independent estimates from the same cohort.

**Fix:** Use ancestry-matched LD scores from alkesgroup (eas_w_ld_chr, afr_w_ld_chr) or compute custom LD scores from in-sample LD reference. For trans-ancestry rg, switch to Popcorn / Popcorn-2.

### Global rg masks local rg variation

**Trigger:** Two traits with biologically plausible shared etiology return global rg ~ 0 in cross-trait LDSC.

**Mechanism:** Global rg averages over the genome; loci with positive local rg can cancel against loci with negative local rg, particularly for traits with antagonistic pleiotropy (e.g. autoimmune vs infectious-disease susceptibility), or when shared etiology is confined to a small fraction of the genome.

**Symptom:** Well-powered (mean chi-square >> 1.02) global LDSC rg near zero with wide CI overlapping zero, while domain biology, prior co-occurrence studies, or shared-pathway analyses strongly suggest shared etiology.

**Fix:** Run LAVA over the standard ~2495 LDetect loci; per-locus Bonferroni-significant rg at any locus is evidence of localized shared etiology. Annotate hit loci with overlapping GWAS catalog signals and pathway/tissue enrichment.

### Cross-population rg below 1 even at causal level

**Trigger:** Computing rg between same-trait GWAS in two ancestries (e.g. EUR T2D vs EAS T2D).

**Mechanism:** Even when the trait has the same biological definition, causal variant identity and effect sizes can differ across populations due to gene-environment interaction, allele-frequency divergence, and population-specific epistasis. Popcorn 2016 demonstrated that rg(cross-pop) < 1 is common and biologically real, not a methodological artifact.

**Symptom:** Trans-ancestry rg point estimate around 0.6-0.9 with CI excluding 1 for a trait expected to be "the same disease".

**Fix:** Use Popcorn or Popcorn-2 (within-pop LDSC is invalid for cross-pop rg); interpret rg(cross-pop) < 1 as quantifying population-specific architecture rather than as bias; report rho_ge (effect correlation) and rho_gi (impact correlation) separately.

### Same-Trait Cross-Cohort rg as Consistency Check

**Use case:** Meta-analysis QC -- two same-trait GWAS (e.g. UKB IBD vs FinnGen IBD) yield rg < 1 with CI excluding 1. This is distinct from the cross-population analog above (handled by Popcorn); here both cohorts are same-ancestry but different studies.

**Interpretation:** (a) population-substructure differences, (b) phenotype-definition heterogeneity (e.g. different ICD coding, self-report vs registry), (c) genuine biology (founder effects in isolates like FinnGen).

**Decision rule:** rg ~ 0.9-1.0 with CI overlapping 1 -> consistent enough to meta-analyze; rg ~ 0.7-0.9 -> moderate heterogeneity, consider sensitivity meta with random effects; rg < 0.7 -> re-examine phenotype definitions before meta-analyzing.

### Low chi-square mean (underpowered GWAS)

**Trigger:** Mean chi-square in either input GWAS is below 1.02 (heuristic LDSC threshold).

**Mechanism:** LDSC, HDL, and LAVA all depend on bivariate moments of Z-scores against LD score; weak signal means the slope is dominated by noise.

**Symptom:** LDSC rg SE > 0.2; intercept estimates fluctuate across chromosome; HDL convergence warnings; LAVA returns p > 0.05 at most loci.

**Fix:** Meta-analyze contributing cohorts to push effective N to ~50k or higher before running rg; if meta-analysis is not feasible, report rg as exploratory; do NOT switch methods to "salvage" power, the problem is upstream of method choice.

### LAVA univariate filter ignored

**Trigger:** Reporting LAVA bivariate rg at a locus where the univariate local h2 is non-significant in one or both traits.

**Mechanism:** LAVA's bivariate test is only valid at loci with detectable local heritability in BOTH traits. Without local h2 signal in at least one trait, the bivariate test is unidentified and can return spurious significant rg.

**Symptom:** LAVA bivariate p < 0.05 at loci where univariate h2 p > 0.05 for one trait; rg estimates near +/- 1 (boundary cases).

**Fix:** Filter loci on `univ.p < 0.05 / N_loci` (Bonferroni for ~2495 loci) in BOTH traits before running `run.bivar()`; report only at filtered loci. This is the documented LAVA workflow in Werme 2022 Supplement and the GitHub vignette.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| |rg| > 0.7 | Operational high correlation | Suggests strong shared genetic architecture; near-universal in pairs like MDD-anxiety or LDL-CHD |
| |rg| 0.3-0.7 | Operational moderate correlation | Common in psychiatric / cardiometabolic trait pairs; routine to flag for joint analysis |
| |rg| < 0.3 | Operational low correlation | Globally weak; may still harbor biologically meaningful local rg via LAVA |
| rg SE < 0.05 | Operational reliable estimate | Above this SE, point estimate is uncertain to 1 decimal place |
| HDL sample overlap < 5% | Ning 2020 Nat Genet (Supplement) | Above this, HDL likelihood is misspecified and biased |
| LDSC mean chi-square > 1.02 | LDSC documentation (Bulik-Sullivan 2015 tutorial) | Below this, LD-score regression is severely underpowered for h2 / rg |
| LAVA local p < 0.05 / N_loci | Werme 2022 Nat Genet 54:274 | Bonferroni for ~2495 LDetect loci; standard genome-wide local-rg correction |
| Popcorn rho_ge CI excludes 1 | Brown 2016 AJHG 99:76 | Evidence of population-specific causal architecture |
| LCV gcp != 0 (two-sided p < 0.05) | O'Connor & Price 2018 Nat Genet 50:1728 | Directional evidence of (partial) causation given non-zero rg |
| LDSC ratio < 0.2 | Bulik-Sullivan 2015 Methods | High ratio (intercept / chi-square - 1) indicates population stratification or model misfit |
| MR + rg sensitivity trigger | Operational | |rg| > 0.3 with a significant IVW estimate REQUIRES CHP-aware sensitivity (CAUSE / LHC-MR) |
| Conditional-rg LAVA covariate set | Werme 2022 | Up to 4 conditioning traits per `run.pcor()` call before identification fails |

## Cross-Trait LDSC: Standard Workflow

**Goal:** Estimate global rg from two GWAS summary statistics, robust to any sample overlap.

**Approach:** Munge each sumstats file (column harmonization + filters), supply ancestry-matched LD scores, run `--rg` mode; interpret slope (rg) and intercept (overlap proxy) separately.

```bash
# Step 1: munge each GWAS to LDSC format (harmonize columns, filter on MAF and INFO, restrict to HapMap3)
munge_sumstats.py \
    --sumstats trait1.tsv.gz \
    --N 250000 \
    --merge-alleles w_hm3.snplist \
    --out trait1.munged

munge_sumstats.py \
    --sumstats trait2.tsv.gz \
    --N 180000 \
    --merge-alleles w_hm3.snplist \
    --out trait2.munged

ldsc.py \
    --rg trait1.munged.sumstats.gz,trait2.munged.sumstats.gz \
    --ref-ld-chr eur_w_ld_chr/ \
    --w-ld-chr eur_w_ld_chr/ \
    --out rg_t1_t2

grep -A 11 'Summary of Genetic Correlation Results' rg_t1_t2.log
```

The log block reports rg, SE, p-value, h2 per trait, and `gcov_int` (genetic covariance intercept = phenotypic-correlation overlap proxy). When running rg of one base trait against many others, use comma-separated lists: `--rg base.sumstats.gz,t1.gz,t2.gz,t3.gz`.

## HDL: Lower-Variance rg for Independent Samples

**Goal:** Estimate global rg with ~60% lower variance than LDSC when the two GWAS draw from non-overlapping samples.

**Approach:** Format each GWAS as an HDL data frame; supply the UKB-array SVD eigen reference path; pass `N0` (overlapping sample count, 0 for independent).

```r
# remotes::install_github('zhenin/HDL/HDL')
library(HDL)

gwas1 <- data.frame(
    SNP = trait1$rsid,
    A1 = trait1$effect_allele,
    A2 = trait1$other_allele,
    N = trait1$N,
    Z = trait1$beta / trait1$se,
    b = trait1$beta,
    se = trait1$se
)

gwas2 <- data.frame(
    SNP = trait2$rsid,
    A1 = trait2$effect_allele,
    A2 = trait2$other_allele,
    N = trait2$N,
    Z = trait2$beta / trait2$se,
    b = trait2$beta,
    se = trait2$se
)

res <- HDL.rg(
    gwas1.df = gwas1,
    gwas2.df = gwas2,
    LD.path = 'UKB_array_SVD_eigen99_extraction',
    N0 = 0,  # 0 if truly independent samples; >0 biases the estimate
    output.file = 'hdl_rg.txt'
)

print(res$rg)
print(res$rg.se)
print(res$P)
```

Pre-download the UKB-array SVD reference (`HDL_documentation.html` -> "How to obtain LD reference panel" link); the eigen99 panel retains 99% of eigenvalues, eigen90 is faster but less precise. Do NOT run HDL when sample overlap is unknown or non-trivial; the wrapper does not warn.

## LAVA: Local Genetic Correlation Per Locus

**Goal:** Identify loci where two traits share genetic etiology, including loci hidden by global rg cancellation.

**Approach:** Process inputs once -> filter to loci with detectable univariate local h2 in BOTH traits -> run bivariate per-locus rg; apply Bonferroni for ~2495 loci.

```r
# remotes::install_github('josefin-werme/LAVA')
library(LAVA)

input <- process.input(
    input.info.file = 'input.info.txt',
    sample.overlap.file = 'sample.overlap.txt',
    ref.prefix = '1kg_EUR_chr',
    phenos = c('trait1', 'trait2')
)

loci <- read.loci('blocks_s2500_m25_f1_w200.GRCh38.locfile')

univ_results <- list()
biv_results <- list()
N_loci <- nrow(loci)

for (i in seq_len(N_loci)) {
    locus <- process.locus(loci[i, ], input)
    if (is.null(locus)) next  # no SNPs / no h2 / monomorphic
    univ <- run.univ(locus)
    univ_results[[i]] <- univ
    pass_univ <- all(univ$p < 0.05 / N_loci)  # Bonferroni on both traits
    if (!pass_univ) next
    biv_results[[i]] <- run.bivar(locus)
}

univ_df <- do.call(rbind, univ_results)
biv_df <- do.call(rbind, biv_results)
biv_df$padj <- p.adjust(biv_df$p, method = 'bonferroni', n = N_loci)
sig_loci <- subset(biv_df, padj < 0.05)
```

The standard LDetect partitioning files (~2495 EUR loci, ~1700 EAS, ~2700 AFR) are at the LAVA GitHub. Sample-overlap file (a per-pair phenotypic-correlation matrix; LDSC intercept is the standard proxy) protects LAVA from the same overlap bias that LDSC's intercept absorbs. For partial / conditional local rg, use `run.pcor(locus, target = c('phenoA', 'phenoB'), phenos = c('cond1', 'cond2', ...))` with up to 4 conditioning traits; the canonical multi-predictor regression alternative is `run.multireg()`.

## rho-HESS: Alternative Local rg

**Goal:** Per-locus bivariate rg using the HESS quadratic-form estimator; complementary to LAVA.

**Approach:** Estimate local h2 per trait first; then bivariate cross-trait estimator using LD-projected effect estimates per locus.

```bash
# Local h2 step (per chromosome)
hess.py \
    --local-hsqg trait1.sumstats.gz \
    --chrom 22 \
    --bfile 1kg_EUR_chr22 \
    --partition-file fourier_ls-chr22.bed \
    --out hess_t1_chr22

# Bivariate local rg (per chromosome, after h2 per trait)
hess.py \
    --local-rhog trait1.sumstats.gz,trait2.sumstats.gz \
    --chrom 22 \
    --bfile 1kg_EUR_chr22 \
    --partition-file fourier_ls-chr22.bed \
    --num-shared 0 \
    --out hess_rhog_chr22
```

`--num-shared` is the number of overlapping individuals; set 0 only if truly independent. HESS partition files use the Berisa & Pickrell 2016 LD-block boundaries (`fourier_ls-*.bed`). LAVA has largely superseded HESS for new analyses, but rho-HESS remains in active use for replication / triangulation.

## Popcorn: Trans-Ancestry rg

**Goal:** Quantify shared causal architecture between two ancestries (e.g. EUR T2D vs EAS T2D) under a MAF-LD-aware model.

**Approach:** Compute cross-population LD scores once per ancestry pair; fit rg using sumstats from each population.

```bash
# Step 1: cross-population LD scores (one-time per ancestry pair)
popcorn compute -v 1 \
    --bfile1 1kg_EUR \
    --bfile2 1kg_EAS \
    --SNPs_to_store hapmap3_snps.txt \
    eur_eas_scores.txt

# Step 2: fit cross-population rg
popcorn fit -v 1 \
    --cfile eur_eas_scores.txt \
    --sfile1 t2d_eur.sumstats.txt \
    --sfile2 t2d_eas.sumstats.txt \
    t2d_eur_eas_rg.txt
```

Popcorn reports `rho_ge` (correlation of causal effect sizes) and `rho_gi` (correlation of variant-level impacts, MAF-weighted). When MAFs differ markedly across populations, `rho_ge` and `rho_gi` diverge; both are biologically meaningful and report-worthy. Effective N per population must be > ~5000 for stable estimates.

## Reconciliation Across Methods

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| LDSC rg and HDL rg agree (independent samples) | Both methods converging on true value | Report HDL as primary (lower SE); LDSC as sensitivity |
| LDSC rg substantially below HDL rg, LDSC intercept large | Sample overlap; HDL is biased toward phenotypic correlation | Report LDSC as primary; flag overlap; do not report HDL |
| Global LDSC rg ~ 0 but LAVA shows multiple Bonferroni-significant local rg | Locus-level cancellation in global average | Report both; the biology is "shared at specific loci, divergent overall" |
| LAVA significant but univariate h2 non-significant at hit locus | Spurious bivariate without identified local h2 signal | Filter univariate first; do NOT report bivariate at unidentified loci |
| Popcorn rho_ge << 1 across many trait pairs | Population-specific causal architecture | Real finding; report rho_ge alongside within-pop h2 |
| LDSC ratio > 0.2 | Population stratification or model misfit | Re-check ancestry; consider LD-score reference mismatch; report with caveat |
| LCV gcp ~ 0 with large rg | Genetic correlation without (partial) causation in either direction | Shared confounder hypothesis is preferred; do NOT report as causal |
| LCV gcp > 0 (significant) with large rg | Partial-to-full causation in tested direction | Combine with bidirectional MR + CHP-aware sensitivity; this is supportive but not sufficient |

**Operational rule for publication:** Report LDSC rg + intercept as primary global statistic; report HDL only if overlap is verified < 5%; complement with LAVA local rg whenever global rg is near zero or biology suggests heterogeneity; report LCV gcp when downstream MR is planned; trans-ancestry analyses require Popcorn (not within-population LDSC).

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| LDSC `category not found` after `--rg` | LD score column header mismatch (custom reference) | Inspect M_5_50 and `.l2.ldscore` headers; align with `--ref-ld-chr` prefix |
| LDSC ratio > 1 (negative h2 z-score) | Severe stratification or wrong LD-score ancestry | Switch to ancestry-matched reference; check for population structure |
| HDL convergence warning / NA SE | Reference panel mismatch or extreme overlap | Verify SVD eigen reference path; switch to LDSC when overlap suspected |
| LAVA `Insufficient SNPs at locus` for most loci | LD reference and partition file from different builds | Match GRCh37 vs GRCh38; align LD reference to partition file |
| LAVA bivariate rg = +/- 1 at boundary | Univariate filter not applied; locus is unidentified | Apply `univ$p < 0.05/N_loci` filter to BOTH traits before `run.bivar()` |
| Popcorn complains about MAF format | sumstats EAF column missing or NA | Provide EAF; do not impute from external reference (creates miscalibration) |
| HESS `--num-shared` defaulting to wrong value | Forgot to set explicitly; default 0 is independent | Always set explicitly; if unknown, use LDSC intercept to infer overlap |
| GenomicSEM `LDSC()` returns negative-definite covariance | Numerical instability with many traits | Inspect per-pair LDSC results; drop low-h2 traits; regularize |
| rg point estimate > 1 with CI overlapping 1 | Sampling variance; same-trait pair near identity | Report as "rg not distinguishable from 1"; constrained likelihood at the rg=1 boundary gives different SE -- LRT against H0: rg=1 is more precise than Wald CI |
| Binary-vs-continuous trait pair scale concern | Reviewer asks about liability-vs-observed scale propagation | LDSC rg is scale-invariant -- case-control h2 liability vs observed scale propagates equivalently into rg; no correction needed |

## Required Reporting for rg Analyses

| Component | Required |
|-----------|----------|
| Per-trait h2 + SE + intercept + mean chi-square | Yes |
| Bivariate rg + SE + p | Yes |
| gcov_int (cross-trait intercept) | Yes; non-zero under known overlap is expected, not bias |
| LD reference panel + ancestry | Yes |
| Method used (LDSC / HDL / LAVA / Popcorn) | Yes; rationale per Decision Tree |
| Local rg supplementary (LAVA) | If global rg null but biology suggests sharing |
| Sample-size: Neff per trait | Yes |

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Sample overlap?" | LDSC: gcov_int reported; non-zero under known overlap is expected, NOT bias. HDL: only valid if independent (<5% overlap) |
| "Why LDSC not HDL?" | HDL gives lower variance but is biased > 5% overlap; LDSC is the conservative default |
| "Local vs global rg?" | If global rg modest but biology suggests sharing, LAVA (Werme 2022) reported as supplementary |
| "Cross-ancestry?" | Popcorn / Popcorn-2 for trans-ancestry; rg < 1 in trans is real biology, not noise |
| "Does rg motivate CHP-MR?" | If `|rg| > 0.3`, CAUSE / LHC-MR sensitivity reported (cross-ref pleiotropy-detection) |
| "rg = 1 boundary?" | If CI includes 1, reported as "not distinguishable from rg=1"; constrained LRT alternative provided |

## Tool Installation Notes

```bash
# LDSC Python 3 fork (original bulik/ldsc is Python 2.7 unmaintained since 2019).
# belowlab/ldsc v3.0.1 broke the --h2/--rg/--h2-cts CLI per its README;
# abdenlab/ldsc-python3 (v2.0.0) retains the working CLI. Docker
# `jtb114/ldsc:latest` is the recommended belowlab fallback.
git clone https://github.com/abdenlab/ldsc-python3.git
conda env create -f ldsc-python3/environment.yml -n ldsc
# pre-computed EUR / EAS / AFR LD scores at alkesgroup.broadinstitute.org/LDSCORE

# HESS
git clone https://github.com/huwenboshi/hess.git
# Berisa-Pickrell LDetect partition files bundled in repo

# Popcorn / Popcorn-2
git clone https://github.com/brielin/Popcorn.git
cd Popcorn && python setup.py install
```

```r
# HDL
remotes::install_github('zhenin/HDL/HDL')
# UKB-array SVD eigen reference: HDL GitHub README has download link

# LAVA
remotes::install_github('josefin-werme/LAVA')
# Pre-computed LDetect partitioning at LAVA GitHub (s2500_m25_f1_w200 for GRCh38, others for GRCh37)

# GenomicSEM
remotes::install_github('GenomicSEM/GenomicSEM')
```

## References

- Bulik-Sullivan B et al 2015 Nat Genet 47:1236 (cross-trait LDSC; intercept absorbs sample overlap)
- Bulik-Sullivan B et al 2015 Nat Genet 47:291 (univariate LDSC h2; companion paper)
- Ning Z et al 2020 Nat Genet 52:859 (HDL; high-definition likelihood; ~60% lower variance than LDSC)
- Werme J et al 2022 Nat Genet 54:274 (LAVA; local genetic correlation via per-locus PC projection)
- Shi H et al 2017 AJHG 101:737 (rho-HESS; locus-level bivariate)
- Shi H et al 2016 AJHG 99:139 (HESS univariate; companion)
- Lee SH et al 2012 Bioinformatics 28:2540 (GREML-bivariate)
- Brown BC et al 2016 AJHG 99:76 (Popcorn; trans-ancestry rg)
- Galinsky KJ et al 2019 PLoS Genet (Popcorn-2 / admixture-aware LD scores)
- Grotzinger AD et al 2019 Nat Hum Behav 3:513 (GenomicSEM)
- O'Connor LJ & Price AL 2018 Nat Genet 50:1728 (LCV; gcp parameter)
- Morrison J et al 2020 Nat Genet 52:740 (CAUSE; CHP-aware MR motivated by high rg)
- Berisa T & Pickrell JK 2016 Bioinformatics 32:283 (LDetect LD blocks underpinning LAVA / HESS)
- Speed D et al 2020 Nat Methods 17:1023 (model comparison of heritability frameworks)
- Bulik-Sullivan B 2015 bioRxiv 018283 (LDSC software documentation and intercept interpretation)
- Skrivankova VW et al 2021 JAMA 326:1614 (STROBE-MR; rg reporting in MR context)

## Related Skills

- causal-genomics/mendelian-randomization - Primary causal estimation; |rg| > 0.3 motivates CHP-aware sensitivity
- causal-genomics/pleiotropy-detection - CAUSE, LHC-MR, LCV; CHP-aware MR battery triggered by high rg
- causal-genomics/heritability-partitioning - Partner method; LDSC stack for univariate h2 and partitioned enrichment
- causal-genomics/genomic-sem - GenomicSEM `ldsc()` is the multivariate extension of bivariate LDSC
- causal-genomics/colocalization-analysis - Locus-level shared causal variant; complements LAVA hits
- causal-genomics/fine-mapping - Credible-set construction at LAVA-significant loci
- causal-genomics/mediation-analysis - MVMR for X -> M -> Y after rg motivates causal hypothesis
- population-genetics/association-testing - GWAS summary statistics underlying all rg methods
- clinical-biostatistics/effect-measures - Translate genetic-architecture findings to clinical effect measures
