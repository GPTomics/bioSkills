---
name: bio-population-genetics-association-testing
description: Single-variant common-variant GWAS with plink2 --glm (linear/logistic, Firth) and the linear mixed models GEMMA, BOLT-LMM, SAIGE, regenie (SPA). A GWAS statistic is valid only when genotype is independent of unmodeled phenotype drivers after the chosen covariates and random effects, so the engine follows sample structure and case:control imbalance, not taste: PC covariates absorb continuous ancestry but cannot remove relatedness (a covariance structure needing an LMM), genomic inflation above 1 is mostly true polygenic signal not confounding (the LDSC intercept is the diagnostic), LOCO prevents proximal contamination, and SPA/Firth keep the tail calibrated at extreme imbalance and low MAC. Use when running single-variant GWAS, choosing between a GLM and a mixed model, or controlling stratification, relatedness, and case:control imbalance. For rare-variant aggregation (burden, SKAT, SKAT-O, ACAT) see rare-variant-association; fine-mapping and MR see causal-genomics; PRS see clinical-databases/polygenic-risk.
tool_type: cli
primary_tool: plink2
---

## Version Compatibility

Reference examples tested with: PLINK 2.0 (alpha 6+), SAIGE 1.3+, regenie 3.4+, BOLT-LMM 2.4+, GEMMA 0.98+, numpy 1.26+, pandas 2.2+, scipy 1.12+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Version traps that change results, not just syntax: PLINK 2.0 `--glm firth-fallback` is the DEFAULT for binary traits and writes `.glm.logistic.hybrid` (mixed logistic and Firth rows), not `.glm.logistic`. SAIGE `--LOCO=TRUE` is the recommended default and the step-1 and step-2 sample IDs plus variance-ratio file must match exactly. regenie applies Firth/SPA only when `--firth`/`--spa` are explicitly set in step 2 (not automatic for every variant). BOLT-LMM is calibrated only for quantitative traits at large N (case fraction >= 10%, MAF > 0.1% for binary coding). The single source of truth for versions is this block, not headings.

# Association Testing

**"Run a GWAS on my genotypes"** -> Fit one regression per variant (additive dosage as predictor) under a confounder model chosen to make genotype independent of unmodeled phenotype drivers, then read effect sizes and p-values that are honest only insofar as that model holds.
- CLI: `plink2 --glm firth-fallback hide-covar --covar pcs.eigenvec` for an unrelated sample whose structure is captured by PCs
- CLI: `SAIGE` (SPA) or `regenie --step 1/2 --firth` for relatedness and/or case:control imbalance; `GEMMA -lmm` or `BOLT-LMM --lmm` for related quantitative traits

Scope: single-variant common-variant GWAS (linear/logistic GLM, linear mixed models, SPA, Firth) and the choice between them. Rare-variant GENE-BASED AGGREGATION (burden, SKAT, SKAT-O, ACAT, STAAR) routes to rare-variant-association. Fine-mapping, Mendelian randomization, and TWAS route to causal-genomics/fine-mapping and causal-genomics/mendelian-randomization. Polygenic risk scores route to clinical-databases/polygenic-risk. Imputed dosages enter from phasing-imputation/genotype-imputation.

## The Single Most Important Insight -- a GWAS statistic is valid only if genotype is independent of unmodeled phenotype drivers AFTER the chosen covariates and random effects

1. Every classic GWAS pathology is one violation of that assumption with a specific signature: structure or cryptic relatedness lifts the whole QQ plot; a heritable covariate adjusted as a collider manufactures direction-flipped hits; case:control imbalance with low MAC under a score test makes the significant tail anti-conservative.
2. The corollary that organizes the toolchain: under a real polygenic trait genomic inflation (lambda) is EXPECTED to exceed 1 and is mostly true signal, so genomic control over-corrects and erases discoveries; the LDSC intercept (not lambda) separates confounding from polygenicity (Bulik-Sullivan 2015).
3. The fix differs per distortion, and the wrong fix makes it worse: lambda-correcting a polygenic trait, or adding PCs to a family sample, both silently degrade the result.
4. Choosing the engine is choosing the assumption: a PC-adjusted GLM assumes structure is a low-rank mean shift, an LMM assumes a polygenic covariance, SPA/Firth assume the score/Wald tail needs the cumulant-generating-function correction.

## Tool Taxonomy

| Method | Citation | Mechanism | When |
|--------|----------|-----------|------|
| plink2 `--glm` | Chang 2015 | Per-variant linear/logistic GLM; Firth fallback for separation | Unrelated sample, common variants, structure captured by PCs |
| GEMMA `-lmm` | Zhou & Stephens 2012 | Exact LMM; fits sigma_g^2 once under the null, then per-variant Wald/LRT/score | Relatedness/structure, quantitative trait, small-medium N |
| BOLT-LMM `--lmm` | Loh 2015 | Variational LMM with a non-infinitesimal Bayesian mixture prior; LD-Score calibrated | Quantitative trait, very large N; gains power with large-effect loci |
| SAIGE | Zhou 2018 | Sparse-GRM LMM + saddlepoint approximation (SPA) on the score statistic | Binary trait with relatedness AND case:control imbalance |
| regenie `--step 1/2` | Mbatchou 2021 | Whole-genome ridge (LOCO) predictor in step 1, Firth/SPA test in step 2, no GRM eigendecomposition | Biobank scale, mixed binary+quantitative; one pipeline |

## Decision Tree by Scenario

| Scenario | Use | Why |
|----------|-----|-----|
| Unrelated, common variants, PCs absorb structure (LDSC intercept ~ 1) | plink2 `--glm` + PC covariates | Fast and exact; an LMM is unnecessary when PCs suffice |
| Relatedness, family, or fine-scale structure | any LMM (GEMMA/BOLT/SAIGE/regenie) | PCs cannot remove a covariance structure; the GRM random effect can |
| Quantitative trait, related, small-medium N | GEMMA `-lmm` or GCTA `--mlma-loco` | Exact LMM; LOCO avoids proximal contamination |
| Quantitative trait, biobank N (>5000) | BOLT-LMM `--lmm` | Scales; non-infinitesimal model adds power; calibrated only at large N |
| Binary trait, related AND imbalanced case:control | SAIGE (SPA) or regenie (`--firth`/`--spa`) | SPA/Firth keep the tail calibrated under imbalance and low MAC |
| Binary trait, unbalanced, NO relatedness | plink2 `--glm firth-fallback` (default) | Firth handles separation; no GRM needed |
| Imputed variants | regress on DOSAGES not hard calls | hard-calling discards imputation uncertainty and biases the SE |
| Rare-variant signal at low MAC | rare-variant-association (burden/SKAT/SKAT-O/ACAT) | single-variant tests are powerless at low MAC; aggregate in a region/gene |

## plink2 GLM (unrelated sample, PC-adjusted)

```bash
# Logistic for binary, linear for quantitative is auto-detected from the phenotype coding.
# firth-fallback is the binary-trait default: ordinary logistic, falling back to Firth only on
# non-convergence (separation). Output is .glm.logistic.hybrid (mixed logistic and Firth rows).
plink2 --bfile qc --pheno pheno.txt --glm firth-fallback hide-covar \
    --covar pcs.eigenvec --covar-name PC1-PC10 --out gwas

# Regress on imputed DOSAGES, not hard calls, so imputation uncertainty enters the SE. Reporting
# columns add A1 frequency and the imputation R2 (machr2 is meaningful only on dosage input) for
# downstream QC and effect-allele bookkeeping.
plink2 --pfile imputed --pheno pheno.txt --glm firth-fallback hide-covar cols=+a1freq,+machr2 \
    --covar covars.txt --covar-name PC1-PC10,age,sex --out gwas_dosage
```

PCs must be computed on LD-pruned, MAF-filtered genotypes with long-range-LD regions excluded (see population-structure); too few PCs leave residual stratification, too many absorb real signal. `--glm sex` adds sex as a covariate on chrX (`no-x-sex` suppresses it); split the pseudoautosomal region before any chrX test (see plink-basics).

## Linear mixed model with LOCO (relatedness / structure)

```bash
# GEMMA: -gk builds the GRM (1=centered, 2=standardized), then -lmm fits one variance component.
# The covariate file passed to -c MUST contain an explicit intercept column of 1s (GEMMA does not add one).
gemma -bfile qc -gk 1 -o grm
gemma -bfile qc -k output/grm.cXX.txt -c covars_with_intercept.txt -lmm 4 -o lmm   # -lmm 4 = Wald+LRT+score

# BOLT-LMM: --lmm decides by cross-validation whether the non-infinitesimal model adds power.
# --lmmInfOnly forces the standard infinitesimal model; --lmmForceNonInf forces the mixture.
bolt --bfile=qc --phenoFile=pheno.txt --phenoCol=trait \
    --covarFile=covars.txt --qCovarCol=PC{1:10} --lmm --LDscoresFile=LDSCORE.1000G_EUR.tab.gz \
    --statsFile=bolt.stats

# SAIGE: step 1 fits the null GLMM (sparse GRM, variance ratio); step 2 runs the SPA score test per variant.
# --LOCO=TRUE (default, recommended) excludes the tested chromosome from the polygenic predictor.
step1_fitNULLGLMM.R --plinkFile=qc --phenoFile=pheno.txt --phenoCol=trait --traitType=binary \
    --covarColList=PC1,PC2,age,sex --sampleIDColinphenoFile=IID --outputPrefix=step1 --LOCO=TRUE
step2_SPAtests.R --vcfFile=chr1.vcf.gz --GMMATmodelFile=step1.rda --varianceRatioFile=step1.varianceRatio.txt \
    --minMAC=20 --is_Firth_beta=TRUE --LOCO=TRUE --SAIGEOutputFile=chr1.saige

# regenie: step 1 builds the whole-genome LOCO ridge predictor; step 2 tests with Firth (--approx for speed) or SPA.
regenie --step 1 --bed qc --phenoFile pheno.txt --covarFile covars.txt --bt --bsize 1000 --out step1
regenie --step 2 --bed qc --phenoFile pheno.txt --covarFile covars.txt --bt \
    --firth --approx --pThresh 0.01 --pred step1_pred.list --bsize 400 --out step2
```

LOCO (leave-one-chromosome-out) is not optional: if the tested variant's chromosome is in the GRM/predictor, the random effect explains part of the variant's own signal and deflates power (proximal contamination). A hand-rolled "GRM from all SNPs" silently throws away power at every true locus.

## Diagnose inflation: lambda vs the LDSC intercept

```python
import numpy as np
from scipy import stats

def lambda_gc(pvalues):
    chisq = stats.chi2.ppf(1 - pvalues, 1)
    return np.median(chisq) / stats.chi2.ppf(0.5, 1)

# lambda > 1 under a polygenic trait is EXPECTED and mostly true signal. Do NOT divide statistics by it.
# Rescale to lambda_1000 (per 1000 cases/1000 controls) before comparing studies of different N.
def lambda_1000(lam, n_cases, n_controls):
    return 1 + (lam - 1) * (1 / n_cases + 1 / n_controls) / (1 / 1000 + 1 / 1000)

# The confounding diagnostic is the LDSC INTERCEPT (run ldsc on the sumstats), not lambda:
# intercept ~ 1 with high lambda = polygenicity (clean); intercept materially > 1 = confounding.
# Prefer the attenuation ratio = (intercept - 1) / (mean(chi2) - 1): the fraction of inflation NOT due
# to polygenicity (~0-0.2 acceptable). The intercept is also inflated by sample overlap in bivariate LDSC.
```

## Per-Method Failure Modes

### Genomic control on a polygenic trait
**Trigger:** dividing every chi-square by lambda_GC because lambda > 1.1. **Mechanism:** lambda rises with N and h2 under true polygenicity, so it is mostly signal. **Symptom:** discoveries vanish; power destroyed. **Fix:** never lambda-correct on lambda alone; use the LDSC intercept and only deflate if intercept-minus-1 is materially > 0.

### PCs for a related sample
**Trigger:** adding "10 PCs" to a family or cryptically-related cohort. **Mechanism:** relatedness is a pairwise covariance, not a low-rank mean shift, so no finite PC set removes it. **Symptom:** inflated, miscalibrated tail statistics despite the PCs. **Fix:** switch to an LMM (GEMMA/BOLT/SAIGE/regenie); more PCs is the wrong lever.

### LMM without LOCO
**Trigger:** building the GRM from all chromosomes including the candidate. **Mechanism:** the variant partly explains itself as a random effect. **Symptom:** deflated statistics, lost power at true loci. **Fix:** use the LOCO variant (`--mlma-loco`, SAIGE `--LOCO=TRUE`, regenie step-1 predictor).

### Wald collapse / score anti-conservatism at imbalance
**Trigger:** plain logistic for a rare or near-monomorphic-in-one-arm variant, or a score test at extreme case:control. **Mechanism:** the MLE diverges (Wald BETA/SE -> 0) or the score null is wrong in the tail. **Symptom:** real rare-variant signal looks non-significant, or the significant tail fills with artifacts. **Fix:** Firth (plink2 firth-fallback, regenie `--firth`) or SPA (SAIGE, regenie `--spa`).

### HWE filtered in cases
**Trigger:** HWE filter on the case-only or combined sample as a discovery QC step. **Mechanism:** a true non-additive disease variant legitimately deviates from HWE in cases. **Symptom:** real associations removed before testing. **Fix:** compute HWE in CONTROLS (or founders) only (see plink-basics).

### Adjusting for a heritable covariate (collider)
**Trigger:** GWAS of a trait while adjusting for a genetically influenced covariate (BMI, smoking). **Mechanism:** opens a non-causal genotype-phenotype path. **Symptom:** spurious, often direction-flipped hits at loci affecting the covariate that replicate within the conditioned design. **Fix:** only adjust for covariates not affected by genotype, or interpret as a different (interventional) question (Aschard 2015).

### BOLT-LMM on a binary trait
**Trigger:** treating BOLT as a logistic engine. **Mechanism:** it runs linear regression on case/control coding. **Symptom:** miscalibration and rare-variant false positives under imbalance. **Fix:** use it only for quantitative traits at case fraction >= 10%, MAF > 0.1%, large N; otherwise SAIGE/regenie.

## Quantitative Thresholds

| Quantity | Typical value | Rationale |
|----------|---------------|-----------|
| Genome-wide significance | p < 5e-8 | Bonferroni over ~1e6 independent common-variant tests in European HapMap LD (Pe'er 2008); ancestry/array-dependent |
| ...for African ancestry | ~1e-8 to 3e-8 | less LD = more independent tests; 5e-8 is too lax |
| ...for WGS / rare-variant scans | ~5e-9 or stricter | far larger effective test count; 5e-8 too permissive |
| Suggestive | p < 1e-5 | follow-up convention, not a calibrated threshold |
| lambda_GC | ~1.0-1.05 fine; 1.05-1.10 inspect; >1.10 investigate | scales with N and h2; pair with lambda_1000 and the LDSC intercept |
| MAC floor (single-variant) | MAC >= 20 (>= 10 with SPA/Firth) | below this even SPA/Firth are unstable; aggregate instead |
| SPA/Firth trigger | case:control more extreme than ~1:10, or low MAC | the score/Wald tail breaks exactly there (SAIGE motivated by ~1:600, Zhou 2018) |
| HWE filter (controls only) | p < 1e-6 | catches genotyping artifacts; never case-only as a discovery filter |

Thresholds are conventions, not laws; inspect distributions and verify current best practice before applying numbers blindly.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `.glm.logistic.hybrid` expected `.glm.logistic` | firth-fallback is the binary default and mixes Firth rows | read the `.hybrid` file; the `FIRTH?` column flags Firth rows |
| 0/1 phenotype gives a null GWAS | PLINK reads 1=control, 2=case by default | pass `--1` for 0/1 coding (see plink-basics) |
| Inflated tail despite 10 PCs | relatedness in the sample | switch to an LMM; PCs cannot remove relatedness |
| Lost power at top loci in an LMM | GRM includes the candidate chromosome | enable LOCO |
| Rare-variant hits look non-significant | Wald collapse under separation | use Firth or SPA |
| Flipped BETA cancels in meta-analysis | effect-allele/strand not harmonized | carry CHR, POS, EA, OA, EAF; resolve A/T and C/G palindromes by frequency or drop |
| Discovery effect size too large downstream | winner's curse | use out-of-sample or shrinkage-corrected effects for PRS/power/MR |
| chrX mis-coded | males hemizygous, PAR diploid | `--glm sex`, split PAR first, handle X-inactivation coding explicitly |
| GEMMA model wrong, no error | `-c` does not add an intercept | the covariate file must contain a column of 1s |
| Fixed-effect pooled OR with I^2 = 80% | heterogeneous true effects | check Cochran's Q / I^2; use random-effects or MR-MEGA for trans-ancestry |

## References

1. Devlin B, Roeder K. Genomic control for association studies. Biometrics 1999; 55(4):997-1004. DOI:10.1111/j.0006-341X.1999.00997.x.
2. Purcell S, Neale B, Todd-Brown K, et al. PLINK: a tool set for whole-genome association and population-based linkage analyses. American Journal of Human Genetics 2007; 81(3):559-575. DOI:10.1086/519795.
3. Pe'er I, Yelensky R, Altshuler D, Daly MJ. Estimation of the multiple testing burden for genomewide association studies of nearly all common variants. Genetic Epidemiology 2008; 32(4):381-385. DOI:10.1002/gepi.20303.
4. Zhou X, Stephens M. Genome-wide efficient mixed-model analysis for association studies. Nature Genetics 2012; 44(7):821-824. DOI:10.1038/ng.2310.
5. Yang J, Lee SH, Goddard ME, Visscher PM. GCTA: a tool for genome-wide complex trait analysis. American Journal of Human Genetics 2011; 88(1):76-82. DOI:10.1016/j.ajhg.2010.11.011.
6. Aschard H, Vilhjalmsson BJ, Joshi AD, Price AL, Kraft P. Adjusting for heritable covariates can bias effect estimates in genome-wide association studies. American Journal of Human Genetics 2015; 96(2):329-339. DOI:10.1016/j.ajhg.2014.12.021.
7. Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ. Second-generation PLINK: rising to the challenge of larger and richer datasets. GigaScience 2015; 4:7. DOI:10.1186/s13742-015-0047-8.
8. Loh PR, Tucker G, Bulik-Sullivan BK, et al. Efficient Bayesian mixed-model analysis increases association power in large cohorts. Nature Genetics 2015; 47(3):284-290. DOI:10.1038/ng.3190.
9. Bulik-Sullivan BK, Loh PR, Finucane HK, et al. LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. Nature Genetics 2015; 47(3):291-295. DOI:10.1038/ng.3211.
10. Zhou W, Nielsen JB, Fritsche LG, et al. Efficiently controlling for case-control imbalance and sample relatedness in large-scale genetic association studies. Nature Genetics 2018; 50(9):1335-1341. DOI:10.1038/s41588-018-0184-y.
11. Mbatchou J, Barnard L, Backman J, et al. Computationally efficient whole-genome regression for quantitative and binary traits. Nature Genetics 2021; 53(7):1097-1103. DOI:10.1038/s41588-021-00870-7.

## Related Skills

- plink-basics - QC, phenotype encoding, and the fileset that enters association
- population-structure - PCA covariates for stratification control
- linkage-disequilibrium - LD pruning before PCA and clumping of GWAS hits
- rare-variant-association - gene-based aggregation (burden, SKAT, SKAT-O, ACAT) below the single-variant MAC floor
- causal-genomics/fine-mapping - from an associated locus to a credible set of causal variants
- causal-genomics/mendelian-randomization - GWAS variants as instruments for causal inference
- clinical-databases/polygenic-risk - PRS built from association sumstats
- phasing-imputation/genotype-imputation - imputed dosages that enter --glm
