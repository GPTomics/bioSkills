---
name: bio-clinical-databases-polygenic-risk
description: Constructs and validates polygenic risk scores using LDpred2-auto, SBayesRC, MegaPRS, PRS-CS, PROSPER, MUSSEL, BridgePRS, JointPRS, PRSmix, or PGS Catalog Calculator with ancestry-aware reference panels (HapMap3, UKB-LD), Pejaver-style calibration, and PRS-RS reporting standards. Use when computing PRS for cohorts, applying Whiffin-style absolute-risk transformation, assessing cross-ancestry portability (Martin 2017 / Ding 2023 continuous ancestry), or auditing PRS manuscripts against the 11-item reviewer checklist.
tool_type: mixed
primary_tool: PGS Catalog Calculator
---

## Version Compatibility

Reference examples tested with: bigsnpr 1.12+ (LDpred2; Privé 2020), PRSice-2 2.3.5+, PRS-CS 1.0.0+ (Ge 2019), gctb 2.5+ (SBayesR/SBayesS/SBayesRC; Zheng 2024), LDAK 6.0+ (MegaPRS; Zhang 2021), pgsc_calc 2.0+ (nf-core; Lambert 2024), Hail 0.2.130+, numpy 1.26+, pandas 2.2+. No general FDA PRS guidance document exists as of May 2026; the operative regulatory text is the August 2025 Federal Register notice on Cancer Predisposition Risk Assessment Systems (Class II device with special controls).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name`
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying. The LDpred2-auto `snp_ldpred2_auto()` signature changed in bigsnpr 1.11+; pin `allow_jump_sign = FALSE` and `shrink_corr = 0.95` explicitly.

# Polygenic Risk Scores; Construction, Calibration, Reporting

**'Compute a PRS for my cohort using these GWAS summary statistics'** -> Match variants to target genotypes, choose method by data availability + trait architecture, derive LD-aware effect estimates, score, normalize by ancestry, transform to absolute risk.

- CLI (recommended one-stop): `pgsc_calc --target target.vcf --pgs_id PGS000001` (nf-core, Lambert 2024)
- R (SOTA single-ancestry): `bigsnpr::snp_ldpred2_auto()` (Privé 2020)
- CLI (multi-ancestry SOTA): `PRS-CSx`, `PROSPER`, `MUSSEL`, `BridgePRS`, `JointPRS`
- CLI (SBayesRC with functional annotations): `gctb --sbayes-rc --bfile target --gwas-summary sumstats.ma`
- CLI (legacy baseline): `PRSice_linux` (clumping + thresholding; still cited for some clinical scores)

## Method Landscape: 2026 Operational Ranking

| Method | Approach | Best for | Fails when |
|--------|----------|----------|-----------|
| **SBayesRC** (Zheng 2024 *Nat Genet*) | Bayesian + 96 functional annotations | EUR; sparse traits | Sumstats LD-incoherent with reference; chain divergence (run --impute-summary first) |
| **MegaPRS** (Zhang 2021 *Nat Commun*) | BLD-LDAK heritability model | EUR; sparser traits | Lacks GCTA-model assumptions; legacy GCTA pipelines |
| **LDpred2-auto** (Privé 2020) | Bayesian + auto-tuning | Polygenic EUR | LD ref mismatch (s > 0.05); allow_jump_sign default True (must pin FALSE) |
| **PRS-CS-auto** (Ge 2019 *Nat Commun*) | Continuous-shrinkage prior | Polygenic EUR | Sparse-trait architecture; HapMap3-restricted variants only |
| **lassosum2** (Privé 2022) | Penalized regression | EUR alternative | Highly polygenic (Bayesian methods better); requires tuning data |
| **C+T** (PRSice-2; Choi 2019) | Clumping + thresholding | Legacy clinical scores (PRS313) | Highly polygenic; Bayesian methods dominate |
| **PROSPER** (Zhang 2024 *Nat Commun*) | Ensemble penalized regression | Multi-ancestry, AFR + others | Single-ancestry; tuning set < 1000 |
| **MUSSEL** (Jin 2024 *Cell Genomics*) | Spike-slab + super-learner | Multi-ancestry; admixed AFR | Single-ancestry; lacks tuning data |
| **JointPRS** (Hu 2025 *Nat Commun*) | Data-adaptive Bayesian | Multi-ancestry; sumstats only | Single-ancestry; very small target |
| **PRS-CSx** (Ruan 2022 *Nat Genet*) | PRS-CS multi-ancestry extension | Multi-ancestry with EUR + non-EUR sumstats | Low causal-variant overlap across ancestries |
| **BridgePRS** (Hoggart 2024 *Nat Genet*) | Ridge-bridge sharing | Low-h^2 AFR / low causal overlap | Standard scenarios (PROSPER/MUSSEL win) |
| **PolyPred / PolyPred+** (Weissbrod 2022) | BOLT-LMM + PolyFun-SuSIE | Multi-ancestry; biobank-scale | Small individual-level data; expensive |

**Citation traps caught by senior PIs:**
- **PRSmix** (Truong 2024); *Cell Genomics*, NOT *Nat Genet*.
- **MUSSEL** (Jin 2024); *Cell Genomics*.
- **PROSPER** (Zhang 2024); *Nat Commun*, NOT *Nat Genet*.
- **Hingorani 2023**; *BMJ Medicine* (NOT main *BMJ*).
- **Mavaddat 2023 BOADICEA update**; *Cancer Epi Biomark Prev*, NOT *Nat Genet*.
- **Mullins 2021** is bipolar disorder, NOT MDD (Howard 2019 *Nat Neurosci* = MDD; Mullins 2021 = BD).

## Multi-Ancestry: The Big Problem

Martin 2017 *Nat Genet* established the 4.5x R^2 attenuation between EUR and AFR. Updates:
- **Martin 2019** *Nat Genet*: "Clinical use of current polygenic risk scores may exacerbate health disparities".
- **Mostafavi 2020** *eLife*: PRS accuracy varies even within a single ancestry due to age, sex, SES, GxE.
- **Ding 2023** *Nature* 618: PGS accuracy decays *continuously* along genetic-ancestry continuum (Pearson r = -0.95 vs PC distance from training data). **The 2026 standard is to report PRS performance vs continuous PC distance, NOT discrete ancestry boxes.**
- **Hou 2023** *Nat Commun*: theoretical decay bounds in MAF/h^2/effect-correlation across ancestries.

Multi-ancestry method choice:
1. Individual-level non-EUR training data + tuning set >= 1000: **PROSPER** or **MUSSEL** (top performance).
2. Summary statistics only + small tuning: **JointPRS** or **PRS-CSx**.
3. Low h^2 / very polygenic / low causal overlap: **BridgePRS**.
4. Functional annotations critical + EUR-dominant: **SBayesRC** (cross-ancestry via `--ldm-eigen`).

## Calibration: The Hingorani Reframing

Khera 2018 *Nat Genet* established the clinical-PRS narrative; **Hingorani 2023** *BMJ Med* is the operative critique:

- HR per SD typically lies between 1.3 and 1.7; **similar to family history alone**.
- Among individuals who develop disease, only ~11% are detected at conventional high-risk PRS threshold; 5% false-positive rate.
- CAD top-2.5% PRS captures 7% of cases; breast-cancer top-quintile captures 6%.
- Wald-Hingorani detection-rate / false-positive-rate ratios approach 10:1 for screening; current PRS achieve 2-3:1.

**Calibration mechanics:**
- Cross-ancestry calibration breaks for the *variance* of the PRS distribution, not just the mean (Sun 2021 *Genome Med*).
- Recalibration: subtract conditional mean given first 4-10 PCs; divide by conditional SD; convert to percentile. **Compute PCs in the test cohort, NOT discovery-cohort PCs**.
- For absolute risk: integrate over external incidence curve (BOADICEA v6 / CanRisk for breast cancer; FOS for CAD).

## PRS-RS Reporting Standards (Wand 2021 *Nature*)

22-item checklist. Reviewer-priority items: cohort independence between development + evaluation (item 13); confounder adjustment in evaluation (item 16); absolute-risk reporting (item 19); ancestry composition of validation (item 21). *Nature Genetics*-tier manuscripts without PRS-RS adherence are rejected at review.

## Decision Tree by Scenario

| Scenario | Recommended path | Why |
|----------|------------------|-----|
| EUR cohort + individual-level data | LDpred2-auto or SBayesRC | SBayesRC integrates functional annotations |
| EUR cohort + sumstats only | MegaPRS or SBayesRC | LDpred2 also viable |
| Highly polygenic trait (height, BMI, education) | LDpred2-auto, PRS-CS-auto | Continuous-shrinkage priors well-suited |
| Sparse trait (lipids, AMD) | MegaPRS, SBayesRC | BLD-LDAK or SBayesR mixture priors |
| Multi-ancestry, large tuning set | PROSPER or MUSSEL | Top performance per benchmarks |
| Multi-ancestry, sumstats + small tuning | JointPRS or PRS-CSx | Joint Bayesian framework |
| Multi-ancestry, target = AFR | MUSSEL or BridgePRS | Best non-EUR performance |
| Combining multiple PGS Catalog scores | PRSmix (single trait) or PRSmix+ (cross-trait) | Elastic-net combination |
| Production score for biobank | `pgsc_calc` Nextflow nf-core | Handles liftover, ancestry, normalization automatically |
| No tuning data available | LDpred2-auto, PRS-CS-auto, JointPRS-auto | Bayesian auto-tuning |
| Clinical reporting | Apply Hingorani-aware framing: absolute-risk transform via external incidence curve | HR per SD = 1.3-1.7 alone is inadequate for screening |

## Standard Workflow: LDpred2-auto

**Goal:** Compute LDpred2-auto PRS from sumstats + target genotypes with appropriate LD reference and sample-overlap detection.

**Approach:** Use bigsnpr R package; match variants strand-aware; compute LD or use prebuilt UK Biobank LD; run `snp_ldpred2_auto()` with shrinkage and jump-sign guards; assess `s` parameter for LD mismatch.

```r
library(bigsnpr)
library(data.table)

# Load target genotypes (PLINK .bed/.bim/.fam -> bigSNP object .rds)
# snp_readBed('target.bed', 'target.rds')
obj <- snp_attach('target.rds')
G <- obj$genotypes
map <- obj$map

# Load GWAS summary stats (standardized format)
sumstats <- fread('gwas_sumstats.txt')
# Required columns: chr, pos, a0 (ref), a1 (effect), beta, beta_se, n_eff, p

# Match variants strand-aware (snp_match handles A/T C/G ambiguity)
df_beta <- snp_match(sumstats, map, strand_flip = TRUE)

# Compute LD correlation matrix (in-sample) OR use prebuilt UKB LD
corr <- snp_cor(G, ind.col = df_beta[['_NUM_ID_']],
                size = 3 / 1000,  # 3 cM window
                ncores = 8)

# LDSC heritability estimate + LD mismatch (s) diagnostic
ldsc_res <- snp_ldsc2(corr, df_beta)
h2_est <- ldsc_res[['h2']]
ldsc_s <- ldsc_res[['int']]  # intercept; large => sample overlap

# LDpred2-auto; run multiple chains in parallel
multi_auto <- snp_ldpred2_auto(
    corr, df_beta,
    h2_init = h2_est,
    vec_p_init = seq_log(1e-4, 0.2, 30),
    burn_in = 500, num_iter = 200,
    allow_jump_sign = FALSE,        # CRITICAL: pin to FALSE to avoid sign artifacts
    shrink_corr = 0.95,             # LD-shrinkage guard
    ncores = 8
)

# Filter divergent chains
beta_auto <- sapply(multi_auto, function(x) x$beta_est)
range_auto <- sapply(multi_auto, function(x) diff(range(x$corr_est)))
keep <- range_auto > (0.95 * quantile(range_auto, 0.95))
beta_final <- rowMeans(beta_auto[, keep, drop = FALSE])

# Score
pred <- big_prodMat(G, beta_final, ind.col = df_beta[['_NUM_ID_']])
```

## SBayesRC Workflow (Functional Annotations)

**Goal:** Compute PRS integrating ~96 functional annotations (baseline-LD v2.2) for +14% R^2 over SBayesR.

**Approach:** GCTB CLI with the SBayesRC algorithm; preprocess sumstats with `--impute-summary`.

```bash
# 1. Impute missing variants in sumstats against reference panel
gctb \
    --sbayes-rc \
    --impute-summary \
    --gwas-summary gwas_sumstats.ma \
    --ldm-eigen ukb_eigen_ld.eigen \
    --annot baselineLD_v2.2 \
    --out sbayesrc_preprocessed.ma

# 2. Run SBayesRC main step
gctb \
    --sbayes-rc \
    --gwas-summary sbayesrc_preprocessed.ma \
    --ldm-eigen ukb_eigen_ld.eigen \
    --annot baselineLD_v2.2 \
    --num-chains 4 --chain-length 25000 --burn-in 5000 \
    --out sbayesrc_run

# 3. Score target genotypes
plink2 --bfile target \
       --score sbayesrc_run.snpRes 2 5 8 header \
       --out sbayesrc_scores
```

## PGS Catalog Calculator (Production Pipeline)

**Goal:** Compute multiple published PGS scores in one pipeline with automatic liftover, ancestry projection, and normalization.

**Approach:** Nextflow nf-core `pgsc_calc` workflow handles GRCh37/38 liftover, PC-projection-based ancestry assignment, ambiguous-SNP filtering, and mean/variance normalization.

```bash
# Calculate multiple PGS for cohort
nextflow run pgscatalog/pgsc_calc \
    --input samplesheet.csv \
    --target_build GRCh38 \
    --pgs_id PGS000004,PGS000019,PGS001775 \
    --run_ancestry resources/pgsc_HGDP+1kGP_v1.tar.zst \
    --outdir pgsc_output/ \
    -profile docker
```

## Multi-Ancestry PRS-CSx

**Goal:** Compute multi-ancestry PRS using ancestry-specific GWAS sumstats jointly via continuous-shrinkage Bayesian framework.

**Approach:** PRS-CSx jointly models multiple ancestries with shared shrinkage prior; outperforms per-ancestry PRS-CS for non-EUR.

```bash
python PRScsx.py \
    --ref_dir=ldblk_1kg \
    --bim_prefix=target \
    --sst_file=eur_sumstats.txt,afr_sumstats.txt,eas_sumstats.txt \
    --n_gwas=200000,30000,50000 \
    --pop=EUR,AFR,EAS \
    --out_dir=prscsx_out --out_name=cohort \
    --phi=1e-2  # tune by trait architecture: 1e-2 polygenic, 1e-4 sparse

# Score each ancestry-specific posterior, then combine on tuning set
for pop in EUR AFR EAS; do
    plink2 --bfile target \
           --score prscsx_out/cohort_${pop}_pst_eff_a1_b0.5_phi1e-02.txt 2 4 6 \
           --out scores_${pop}
done
```

## Score Normalization and Ancestry Recalibration

**Goal:** Convert raw PRS to ancestry-conditional percentiles for clinical interpretation.

**Approach:** Subtract conditional mean given test-cohort PCs (NOT discovery PCs); divide by conditional SD.

```python
import numpy as np
from scipy import stats
import statsmodels.api as sm

def ancestry_conditional_normalize(prs, pcs, n_pcs=10):
    '''Recalibrate PRS by removing ancestry effects (Sun 2021 + Ding 2023 continuous).

    pcs: matrix of principal components computed in the TEST cohort (not discovery).
    Returns: ancestry-conditional Z scores.
    '''
    X = sm.add_constant(pcs[:, :n_pcs])
    mean_model = sm.OLS(prs, X).fit()
    expected = mean_model.predict(X)
    residuals = prs - expected

    log_var_model = sm.OLS(np.log(residuals ** 2 + 1e-12), X).fit()
    expected_log_var = log_var_model.predict(X)
    sd = np.sqrt(np.exp(expected_log_var))

    return residuals / sd


def prs_to_percentile(z_scores):
    '''Convert ancestry-conditional Z to population percentile.'''
    return stats.norm.cdf(z_scores) * 100


def absolute_risk_transform(percentile, incidence_curve_age, age):
    '''Integrate over external age-conditional incidence curve to get absolute risk.

    Khera 2018 used FOS for CAD; Lee 2019 BOADICEA v6 for breast cancer.
    The HR-per-SD-only framing (Hingorani 2023 critique) is insufficient for screening.
    '''
    base_risk = incidence_curve_age(age)
    z = stats.norm.ppf(percentile / 100)
    relative_risk = np.exp(z * 0.5)  # placeholder; trait-specific log HR
    return base_risk * relative_risk
```

## Sample Overlap Detection (EraSOR / Bivariate LDSC Intercept)

**Goal:** Detect overlap between PRS-discovery cohort and target cohort, which inflates apparent PRS performance.

**Approach:** Run bivariate LDSC; |intercept| > 0.05 with target n >= 1000 is the typical alarm.

```bash
# bivariate LDSC for sample overlap
ldsc.py \
    --rg target_sumstats.sumstats.gz,discovery_sumstats.sumstats.gz \
    --ref-ld-chr eur_w_ld_chr/ \
    --w-ld-chr eur_w_ld_chr/ \
    --out overlap_check

# Inspect *.log: gcov_int is the sample-overlap-driven intercept (after gencov correction)
# |gcov_int| > 0.05 with target n >= 1000 => substantial overlap
```

## Per-Operation Failure Modes

**1. Discovery + tuning + test overlap (the inflated R^2 trap)**
- Trigger: GWAS run on UKB; PRS evaluated in UKB; tuning on UKB subset.
- Mechanism: Discovery samples leak into test set; PRS-derived effects fit the test data perfectly.
- Symptom: Reported R^2 inflated 2-5x; performance does not replicate in external cohort.
- Fix: Three-way disjoint sample partition (Wray 2014 *J Child Psychol Psychiatry*; Choi 2020 *Nat Protoc*); use **EraSOR** (Choi 2023 *GigaScience*) or bivariate LDSC intercept; threshold |intercept| > 0.05 with target n >= 1000.

**2. Treating discrete ancestry boxes as ground truth**
- Trigger: Report "PRS for AFR cohort" without quantifying genetic distance.
- Mechanism: Ancestry is continuous (Ding 2023 *Nature*); per-PC distance from training data drives R^2 with r = -0.95.
- Symptom: Reviewers reject for incomplete ancestry characterization.
- Fix: Project test cohort onto reference PCs; report performance as a function of continuous PC distance.

**3. Strand-ambiguous SNPs (A/T, C/G) handled wrong**
- Trigger: Drop or trust strand annotation across cohorts.
- Mechanism: ~10-15% of SNPs are strand-ambiguous; cross-platform/cohort scoring requires explicit handling.
- Symptom: Strand-flipped variants reverse effect direction; PRS uncorrelated with phenotype.
- Fix: PRSice-2 default drops them; LDpred2 `snp_match()` frequency-matches with 0.4-0.6 MAF tolerance.

**4. LD reference mismatch (high `s` parameter)**
- Trigger: Use 1KG-EUR LD for UK Biobank target; LDpred2 silently uses mismatched LD.
- Mechanism: `snp_ldsc2()` returns LD-mismatch parameter `s`; values > 0.05 indicate problematic mismatch.
- Symptom: Posterior effect estimates inflated; PRS doesn't replicate.
- Fix: Use UKB LD panel (n=40k+) instead of 1KG-EUR (n=489); inspect `s` routinely.

**5. PCs in derivation but not in evaluation**
- Trigger: GWAS uses 10 PCs; PRS evaluation regression omits PCs.
- Mechanism: Population stratification effects re-enter the PRS-phenotype association without PC control.
- Symptom: Apparent PRS association is partially ancestry confounding, not biology.
- Fix: Include same PCs in evaluation regression; compute PCs in TEST cohort (not discovery).

**6. Treating HR per SD = 1.5 as clinically actionable**
- Trigger: Report "top 5% PRS = 1.5x risk -> screening criterion".
- Mechanism: HR per SD 1.3-1.7 is similar to family history; Wald-Hingorani DR-to-FPR ratios approach 10:1 for true screening utility; PRS achieve 2-3:1.
- Symptom: Clinical implementation produces high false-positive rates.
- Fix: Integrate over external incidence curve for absolute-risk reporting; benchmark against family history specifically.

**7. PRSmix / MUSSEL / PROSPER cited in wrong journal**
- Trigger: Manuscript cites these as *Nat Genet*.
- Mechanism: PRSmix and MUSSEL = *Cell Genomics*; PROSPER = *Nat Commun*.
- Symptom: Senior PIs catch in review.
- Fix: Verify citations against PubMed before submission.

**8. HLA region included in PRS without explicit handling**
- Trigger: Compute autoimmune-disease PRS with HLA region included naively.
- Mechanism: HLA region (chr6 25-35 Mb) has extreme LD; SNP-based PRS does not capture classical HLA allele effects.
- Symptom: Reviewer flags incomplete HLA modeling.
- Fix: Exclude HLA region from main PRS; model classical HLA alleles separately via HIBAG / SNP2HLA / HLA-LA; report HLA-excluded + HLA-augmented versions.

**9. Sex chromosomes mishandled**
- Trigger: PRS includes chrX without sex-specific dosage coding.
- Mechanism: Males 0/1; females 0/1/2 unless XCI-corrected.
- Symptom: Cross-sex comparison invalid.
- Fix: Drop chrX OR apply explicit sex-stratified dosage encoding.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|-------------|--------|
| LDpred2-auto vs PRS-CS-auto disagree on top decile | Trait sparsity differs; method assumption mismatch | Use SBayesRC as tiebreaker; ensemble via PRSmix |
| Multi-ancestry methods rank differently | Tuning set size; ancestry composition | Report ensemble; document hyperparameter sensitivity |
| EUR PRS R^2 << non-EUR target R^2 | Expected; Martin 2017 transferability | Switch to PROSPER / MUSSEL / BridgePRS |
| PRS performance differs across age strata | Mostafavi 2020; within-ancestry heterogeneity | Report age-stratified PRS performance |
| PRS unrelated to phenotype despite high h^2 | Sample overlap; strand-ambiguous SNPs; PC issues | Run EraSOR; check `snp_match()` flips; recompute PCs in test |
| PRS percentile shifts with array platform | Different ascertainment of common variants | Recalibrate per-platform; use ancestry-conditional Z |

## Quantitative Thresholds and Conventions

| Threshold | Convention | Source |
|-----------|-----------|--------|
| HapMap3 SNPs | ~1.1M variants; required for PRS-CS, LDpred2 | HapMap project |
| LDpred2 LD ref `s` | < 0.05 indicates well-matched LD | Privé 2022 |
| EraSOR | |intercept| > 0.05 with n >= 1000 => sample overlap | Choi 2023 |
| MAF range for ambiguous SNPs | Drop or frequency-match 0.4-0.6 | LDpred2 default |
| INFO score (imputed) | >= 0.8 for inclusion | Yengo 2018 |
| Kinship coefficient cutoff | KING > 0.0884 = 3rd-degree or closer; remove | KING documentation |
| HLA region exclusion | chr6 28-34 Mb (some use 25-35) | Convention |
| HR per SD typical | 1.3-1.7 for most diseases | Hingorani 2023 BMJ Med |
| Top 2.5% PRS CAD detection rate | ~7% of cases captured | Hingorani 2023 |
| Mavaddat PRS313 | 313 SNPs, breast cancer | Mavaddat 2019 |
| Khera CAD PRS | ~6.6M variants | Khera 2018 |
| Aragam CAD PRS | GPS_Mult; multi-ancestry SOTA 2023 | Aragam 2023 |

## Common Errors

| Symptom | Cause | Solution |
|---------|-------|----------|
| PRS R^2 << expected on held-out data | Sample overlap with discovery | Run EraSOR; rebuild three-way disjoint splits |
| LDpred2-auto chains diverge | LD mismatch (`s` > 0.05) or wrong h^2 init | Use UKB LD panel; check `snp_ldsc2()` output |
| Strand-flip warnings ignored | `snp_match` defaults | Set `strand_flip = TRUE`; review flipped variants |
| HLA region produces extreme PRS values | Long-range LD | Exclude chr6 28-34 Mb; model HLA separately |
| Non-EUR PRS at chance | Method not multi-ancestry-aware | Switch to PROSPER/MUSSEL/BridgePRS |
| Different ranking with same data | Stochastic MCMC | Seed; report posterior CIs; check convergence |
| PGS Catalog liftover failure | Incompatible build | Use pgsc_calc auto-liftover; check log |
| Top 1% PRS appears benign | Conditional mean drift across ancestries | Apply Sun 2021 / Ding 2023 ancestry-conditional Z |

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "PRS HR 1.5 per SD is similar to family history" | Acknowledged (Hingorani 2023). We report absolute-risk integrated over age-conditional incidence (Wald-Hingorani framing); we do NOT recommend population-level screening on PRS alone. |
| "Why use Bayesian methods over C+T?" | LDpred2-auto / PRS-CS / SBayesRC capture polygenic architecture better; +20-50% R^2 vs C+T per Pain 2021 / Privé 2022 benchmarks. |
| "Why not use the largest published GWAS?" | We checked sample overlap between discovery and target with EraSOR / bivariate LDSC intercept; |intercept| < 0.05. |
| "Discrete ancestry boxes don't capture genetic structure" | Agreed; we report PRS performance vs continuous PC distance (Ding 2023 *Nature*); discrete labels are summary only. |
| "Why exclude HLA?" | HLA region extreme LD makes SNP-based posteriors unstable; we model classical HLA alleles separately via HIBAG / SNP2HLA. |
| "PRS-CS phi parameter; how chosen?" | We use --phi=auto with multiple chains; for sparse traits (lipids) 1e-4; for highly polygenic (height) 1e-2. |
| "Where is the absolute risk?" | Reported per Wand 2021 PRS-RS item 19; integrated over age-conditional incidence using external curve. |
| "FDA PRS guidance?" | No general FDA PRS draft guidance exists as of 2026; August 2025 Federal Register Class II classification for Cancer Predisposition Risk Assessment Systems is the operative regulatory text. |
| "Why three different references for SBayesRC?" | Zheng 2024 *Nat Genet* is the primary paper; baseline-LD v2.2 is Finucane 2018; UKB LD panel is the bigsnpr/Privé 2020 release. |

## References

- Privé F et al. 2020. LDpred2: better, faster, stronger. *Bioinformatics* 36:5424.
- Privé F et al. 2022. Identifying and correcting for misspecifications in GWAS summary statistics and polygenic scores. *Hum Genet Genom Adv* 3:100136. (LDpred2 misspecification + lassosum2)
- Ge T et al. 2019. Polygenic prediction via Bayesian regression and continuous shrinkage priors. *Nat Commun* 10:1776. (PRS-CS)
- Ruan Y et al. 2022. Improving polygenic prediction in ancestrally diverse populations. *Nat Genet* 54:573. (PRS-CSx)
- Zheng Z et al. 2024. Leveraging functional genomic annotations and genome coverage to improve polygenic prediction of complex traits within and between ancestries. *Nat Genet* 56:767. (SBayesRC)
- Lloyd-Jones LR et al. 2019. Improved polygenic prediction by Bayesian multiple regression on summary statistics. *Nat Commun* 10:5086. (SBayesR)
- Zeng J et al. 2021. Widespread signatures of natural selection across human complex traits and major risk factor categories. *Nat Commun* 12:1164. (SBayesS)
- Zhang Q et al. 2021. Improved genetic prediction of complex traits from individual-level data or summary statistics. *Nat Commun* 12:4192. (MegaPRS)
- Zhang H et al. 2024. PROSPER: enhanced polygenic risk score with summary statistics. *Nat Commun* 15:3413.
- Jin Y et al. 2024. MUSSEL: enhanced Bayesian polygenic risk prediction in admixed populations. *Cell Genomics* 4:100539.
- Hu S et al. 2025. JointPRS: a data-adaptive framework for joint polygenic risk score modeling. *Nat Commun* 16:59243.
- Hoggart CJ et al. 2024. BridgePRS leverages shared genetic effects across ancestries. *Nat Genet* 56:180.
- Weissbrod O et al. 2022. Leveraging fine-mapping and multipopulation training data to improve cross-population polygenic risk scores. *Nat Genet* 54:450. (PolyPred)
- Zhao Z et al. 2022. PUMAS: fine-tuning polygenic risk scores with GWAS summary statistics. *Am J Hum Genet* 109:2253. (TL-PRS-adjacent transfer learning)
- Truong B et al. 2024. PRSmix: integrating multiple polygenic risk scores. *Cell Genomics* 4:100495.
- Mostafavi H et al. 2020. Variable prediction accuracy of polygenic scores within an ancestry group. *eLife* 9:e48376.
- Ding Y et al. 2023. Polygenic scoring accuracy varies across the genetic ancestry continuum. *Nature* 618:774.
- Hou K et al. 2023. Causal effects on complex traits are similar for common variants across segments of different continental ancestries within admixed individuals. *Nat Commun* 14:5566.
- Hingorani AD et al. 2023. Performance of polygenic risk scores in screening, prediction, and risk stratification. *BMJ Med* 2:e000554.
- Wand H et al. 2021. Improving reporting standards for polygenic scores in risk prediction studies. *Nature* 591:211. (PRS-RS)
- Lambert SA et al. 2021. The Polygenic Score Catalog as an open database for reproducibility and systematic evaluation. *Nat Genet* 53:420.
- Lambert SA et al. 2024. The PGS Catalog Calculator: automated polygenic scoring of large datasets at scale. *Nat Genet* 56:1989.
- Fritsche LG et al. 2020. Cancer PRSweb: an online repository with polygenic risk scores for major cancer traits and their evaluation in two cohorts. *Am J Hum Genet* 107:815. (PRSweb)
- Khera AV et al. 2018. Genome-wide polygenic scores for common diseases identify individuals with risk equivalent to monogenic mutations. *Nat Genet* 50:1219.
- Aragam KG et al. 2022. Discovery and systematic characterization of risk variants and genes for coronary artery disease in over a million participants. *Nat Genet* 54:1803.
- Aragam KG et al. 2023. Genome-wide polygenic scores reduce risk reclassification in primary prevention. *Nat Med* 29:1793.
- Mavaddat N et al. 2019. Polygenic risk scores for prediction of breast cancer and breast cancer subtypes. *Am J Hum Genet* 104:21. (PRS313)
- Conti DV et al. 2021. Trans-ancestry genome-wide association meta-analysis of prostate cancer identifies new susceptibility loci. *Nat Genet* 53:65.
- Trubetskoy V et al. 2022. Mapping genomic loci implicates genes and synaptic biology in schizophrenia. *Nature* 604:502. (PGC3)
- Howard DM et al. 2019. Genome-wide meta-analysis of depression identifies 102 independent variants. *Nat Neurosci* 22:343. (MDD)
- Mullins N et al. 2021. Genome-wide association study of more than 40,000 bipolar disorder cases. *Nat Genet* 53:817. (BIPOLAR; not MDD)
- Bellenguez C et al. 2022. New insights into the genetic etiology of Alzheimer's disease. *Nat Genet* 54:412.
- Mahajan A et al. 2022. Multi-ancestry genetic study of type 2 diabetes highlights the power of diverse populations. *Nat Genet* 54:560. (DIAMANTE)
- Suzuki K et al. 2024. Genetic drivers of heterogeneity in type 2 diabetes pathophysiology. *Nature* 627:347. (Largest T2D multi-ancestry GWAS)
- PGS Catalog: `https://www.pgscatalog.org`
- PGS Catalog Calculator: `https://github.com/PGScatalog/pgsc_calc`
- Federal Register August 2025 Cancer Predisposition Risk Assessment System: `https://www.federalregister.gov/documents/2025/08/21/2025-16035/`

## Related Skills

- clinical-databases/gnomad-frequencies - Population AF for QC
- clinical-databases/variant-prioritization - Rare-variant filtering background
- clinical-databases/clinvar-lookup - Variant pathogenicity
- causal-genomics/mendelian-randomization - PRS as instrument
- population-genetics/population-structure - Ancestry inference for PC computation
- machine-learning/biomarker-discovery - PRS as biomarker component
