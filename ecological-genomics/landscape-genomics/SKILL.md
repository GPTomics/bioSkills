---
name: bio-ecological-genomics-landscape-genomics
description: Tests genotype-environment associations and identifies adaptive loci while correcting for the four-confound landscape (structure, demography, background selection, sampling design) using LFMM2 with mandatory K via sNMF cross-entropy elbow (LEA 3), BayPass Core/AUX/C2/IS with Omega covariance matrix, RDA / pRDA for polygenic adaptation (Forester 2018; requires imputed genotypes), OutFLANK with trimmed FST null, pcadapt, gradient forests (Ellis-Smith-Pitcher 2012, NOT mis-cited Ellis-Manel), Capblancq & Forester 2021 RDA Swiss-army-knife, genomic-offset prediction with Lind & Lotterhos 2025 three-regime caveat, Lotterhos-Whitlock sampling optima, Wang & Bradburd 2014 IBD vs IBE, and Circuitscape + ResistanceGA. Use when identifying adaptive loci across gradients, choosing K for LFMM2, deciding among GEA methods, predicting maladaptation with the novel-environment caveat, distinguishing IBD vs IBE, or optimizing sampling design.
tool_type: r
primary_tool: LEA
---

## Version Compatibility

Reference examples tested with: LEA 3.14+ (lfmm2 via LEA 3), pcadapt 4.3+, OutFLANK 0.2+, vegan 2.6+, gradientForest 0.1-32+, terra 1.7+, qvalue 2.34+, BayPass 2.3+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Landscape Genomics

**"Find loci associated with environmental adaptation in my populations"** -> Genotype-environment association with K-selected latent-factor correction, multivariate RDA for polygenic signal, paired with demography-calibrated FST outliers; genomic-offset prediction with explicit prediction-novelty regime characterization.
- R: `LEA::lfmm2()` (modern; via LEA 3) for univariate GEA with cross-entropy-elbow K choice
- R: `vegan::rda()` with `Condition()` for partial-RDA polygenic GEA
- R: `OutFLANK::OutFLANK()` for demography-calibrated FST outliers
- R: `pcadapt::pcadapt()` for PC-based scans without environmental data
- CLI: `g_baypass` for BayPass Bayesian Core/AUX/C2/IS analyses

## The Single Most Important Modern Insight -- GEA is a Four-Confound Problem, Not a Signal-Plus-Noise Problem

The number-one failure mode in landscape genomics is treating "is locus X associated with environment Y" as a simple signal-extraction problem (the NimBios working group consolidated this view in Hoban 2016 *Am Nat* 188:379-397). **GEA has FOUR concurrent confounds**:
1. **Population structure** — IBD-driven allele-frequency gradients track distance, which correlates with most environmental variables
2. **Demographic history** — range expansion creates allele-frequency clines along the expansion axis that mimic selection (Lotterhos & Whitlock 2014 *Mol Ecol* 23:2178-2192)
3. **Background and linked selection** — reduces diversity in low-recombination regions, inflating apparent local-FST signals
4. **Sampling design** — non-random spatial sampling creates spurious GEA signals (Lotterhos & Whitlock 2015 *Mol Ecol* 24:1031-1046; **paired contrasts > random > transects** when demography is concerning)

A second cornerstone: **Forester et al. 2018 established RDA as the best polygenic-adaptation detector** but with a caveat — RDA performance requires imputed genotypes (no missing data) and is best for linear-gradient environments. For monogenic strong-selection loci, OutFLANK and BayPass-C2 remain competitive.

A third: **Genomic offset has three regimes** (Lind & Lotterhos 2025 *Mol Ecol Resour* 25:e14008): works well for interpolation (similar-to-training environments), degrades gracefully for modest extrapolation, FAILS for highly novel future environments. Reporting a single offset map without characterizing prediction-novelty is methodologically incomplete. For the broader eco-evolutionary integration framework, see Aguirre-Liguori 2021 *Nat Ecol Evol* 5:1350-1360.

## Algorithmic Taxonomy

| Method | Type | Strength | Fails when |
|--------|------|----------|------------|
| FST outliers (naive) | Univariate per-locus | Simple; widely understood | Inflated FDR under demography (Lotterhos & Whitlock 2014) |
| OutFLANK (Whitlock & Lotterhos 2015) | Univariate FST with trimmed null | Robust to demography via trimmed-tail null | Conservative; misses weak-effect polygenic |
| pcadapt (Luu 2017) | PC-based Mahalanobis | Handles continuous structure; no need to define populations | Detects axes-of-divergence loci; not environment-specific |
| LFMM (Frichot 2013) | MCMC mixed model with latent factors | Original framework | Slow; superseded by LFMM2 |
| LFMM2 (Caye 2019) | LSE-based mixed model with K latent factors | Fast; orders of magnitude faster than LFMM; modern default | K choice is critical; wrong K silently invalidates results |
| BayPass Core (Gautier 2015) | Bayesian FST with covariance matrix Omega | Explicit shared-history correction via Omega | Computationally heavy; convergence needs care |
| BayPass AUX | Bayesian env association with binary auxiliary variable | Posterior gives Bayes Factor per locus | Same as Core |
| BayPass C2 | Bayesian contrast between two pre-defined groups | Like Bayesian Fisher exact | Requires pre-defined groups, not continuous env |
| BayPass IS | Importance-sampling joint Bayesian | Multiple covariates simultaneously | Most computationally expensive |
| RDA (Forester 2018) | Multivariate constrained ordination | HIGH power and LOW FDR for polygenic adaptation | Requires imputed genotypes; less power for monogenic |
| pRDA (partial RDA) | RDA conditioning out population structure | Controls for structure explicitly | Conservative; depends on which axes are conditioned out |
| Gradient forests (Ellis 2012; Ellis-Smith-Pitcher) | Random-forest non-linear GEA | Detects non-linear gene-environment | Less interpretable; no formal p-values |
| RDA offset (Capblancq 2021) | Linear genomic-offset prediction | Linear gradient assumption | Fails under highly novel future climate (Lind 2025) |
| Gradient-forest offset | Non-linear genomic-offset prediction | Captures non-linearity | Same three-regime caveat |
| RONA | Locus-by-locus offset | Simple, transparent | Locus-specific; ignores multilocus structure |
| Circuitscape (McRae 2008) | Circuit-theory landscape resistance | Integrates ALL paths; not single-best | Symmetric only; asymmetric landscapes need directional methods |
| ResistanceGA (Peterman 2018) | Genetic-algorithm optimization of resistance | Optimizes cost surface itself | MLPE parameterization needed (specific lmer form) |

## Decision Tree by Scenario

| Scenario | Recommended approach | Why |
|----------|---------------------|-----|
| Monogenic local adaptation | OutFLANK or BayPass Core; cross-check with simulations | Demography-calibrated null |
| Polygenic local adaptation | RDA with environmental matrix; pRDA conditioning on structure | Forester 2018 superior power + lower FDR |
| Environmental association with structure correction | LFMM2 (K from sNMF cross-entropy elbow) | Modern default; orders of magnitude faster than LFMM |
| Binary group contrast | BayPass C2 | Like Bayesian Fisher exact |
| Joint multi-covariate Bayesian inference | BayPass IS | Multiple env simultaneously |
| Quantifying maladaptation under climate change | RDA offset OR gradientForest offset | Report Lind 2025 caveat: works for interpolation, fails for novel |
| Multi-method consensus | Run OutFLANK + LFMM2 + BayPass; overlap is high-confidence | Method differences = methodological coverage |
| Mapping landscape resistance | Circuitscape + ResistanceGA optimization | Circuitscape integrates all paths; ResistanceGA optimizes the surface |
| Distinguishing IBD vs IBE | dbRDA variance partitioning (geographic + env) | Mantel-based methods have known autocorrelation issues |
| Sampling design optimization | Paired contrasts at environmental endpoints; cite Lotterhos & Whitlock 2015 | Most powerful under demographic concern |
| Polyploid species | polyRAD / fitPoly / updog FIRST; then standard GEA | Diploid-assuming tools give biased FST estimates |

## LFMM2 with Mandatory K Selection

**Goal:** Identify loci significantly associated with environmental variables while controlling for unobserved demographic structure via K latent factors.

**Approach:** Run sNMF on the genotype matrix across K = 1...10 (or higher), compute cross-entropy at each K, pick K at the elbow (where cross-entropy first plateaus), then pass that K to `LEA::lfmm2()`. Wrong K silently invalidates results — too few K means residual structure confounds environment; too many K absorbs the environmental signal.

```r
library(LEA)
library(qvalue)

# Convert VCF -> LFMM/GENO format
vcf2lfmm('variants.vcf', 'genotypes.lfmm')
vcf2geno('variants.vcf', 'genotypes.geno')

# MANDATORY K selection via sNMF cross-entropy elbow
# Run multiple repetitions per K (5 minimum) for stability
snmf_result <- snmf('genotypes.geno', K = 1:10, repetitions = 5,
                     entropy = TRUE, project = 'new')
ce_values <- sapply(1:10, function(k) min(cross.entropy(snmf_result, K = k)))
plot(1:10, ce_values, xlab = 'K', ylab = 'Cross-entropy',
     pch = 19, col = 'blue', type = 'b')

# Pick K at elbow (first plateau); report sensitivity across K
best_K <- which.min(ce_values)

# Run LFMM2 with selected K
genotypes <- read.lfmm('genotypes.lfmm')
env_vars <- read.env('environment.env')
lfmm_result <- lfmm2(input = genotypes, env = env_vars, K = best_K)

# Test associations with genomic-control calibration
pvalues <- lfmm2.test(lfmm_result, input = genotypes, env = env_vars,
                       full = TRUE, genomic.control = TRUE)

# Check genomic inflation factor (GIF / lambda)
# Target ~1.0; > 1.5 means insufficient structure correction (try larger K)
# < 0.5 means over-correction (try smaller K)
gif <- median(qchisq(1 - pvalues$pvalues[, 1], df = 1)) / qchisq(0.5, df = 1)
cat('Genomic inflation factor (lambda):', round(gif, 3), '\n')

# Storey FDR control
qvals <- qvalue(pvalues$pvalues[, 1])$qvalues
candidates <- which(qvals < 0.05)
cat('Candidate adaptive loci (q < 0.05):', length(candidates), '\n')

# Sensitivity check: re-run with K +/- 1 and report overlap of candidates
# Loci detected only at one K are sensitive to latent-factor choice; flag as lower-confidence
```

## RDA — Capblancq & Forester 2021 Swiss-Army-Knife Workflow

**Goal:** Detect polygenic local adaptation via multivariate constrained ordination, with optional climate-offset prediction.

**Approach:** Apply the Capblancq & Forester 2021 RDA workflow: (1) variable selection via forward selection with adjusted R^2; (2) variance partitioning into pure-environment / pure-spatial / shared; (3) GEA outlier identification at SD > 3 on RDA axes 1-3; (4) adaptive-index computation; (5) genomic-offset prediction with the Lind 2025 three-regime caveat. **CRITICAL: genotypes must be imputed before RDA** (no missing data); cite Forester 2018 for the polygenic-power result.

```r
library(vegan)

# CRITICAL: impute missing genotypes BEFORE RDA (Forester 2018 benchmark requirement)
# Column-mean imputation per locus is the simplest valid choice;
# population-stratified mean imputation is preferred when populations are known.
# For LEA-generated data, LEA::impute() is the canonical option.
for (j in seq_len(ncol(allele_freq))) {
    col_mean <- mean(allele_freq[, j], na.rm = TRUE)
    allele_freq[is.na(allele_freq[, j]), j] <- col_mean
}

# Partial RDA conditioning on population structure (e.g., Q-matrix from sNMF)
# This is pRDA per Capblancq 2021
rda_result <- rda(allele_freq ~ temperature + precipitation + altitude +
                  Condition(as.matrix(q_matrix)), data = env_data)

# Permutation test
anova(rda_result, permutations = 999)

# Variance partitioning (Capblancq workflow step 2)
vp <- varpart(allele_freq, env_data, q_matrix)
vp$part$fract  # pure env / shared / pure structure / residual

# GEA outliers (z-score > 3 on RDA axes 1-3)
# Threshold of 3 SD is conservative; for polygenic adaptation, may use 2.5 SD
loadings <- scores(rda_result, choices = 1:3, display = 'species')
zscores <- apply(loadings, 2, function(x) (x - mean(x)) / sd(x))
rda_candidates <- which(apply(abs(zscores), 1, max) > 3)
cat('RDA candidate loci (|z| > 3):', length(rda_candidates), '\n')
```

## BayPass — Core, AUX, C2, IS Models

**Goal:** Bayesian Genotype-Environment Association with explicit population-history correction via covariance matrix Omega.

**Approach:** Convert VCF to BayPass's space-separated allele-count format. Run Core (estimates Omega; like Bayesian FST); use AUX for binary env-association testing; C2 for two-group contrast; IS for joint multi-covariate.

```bash
# BayPass input: space-separated allele counts, NOT VCF
# Use vcf2baypass.pl OR bcftools query to convert

# Core model: estimate Omega (population covariance from genome-wide SNPs)
g_baypass -gfile geno_BayPass.txt -outprefix core_run -nthreads 4

# AUX (binary auxiliary variable): test environmental association
# omegafile from Core run; covariates is samples x env file
g_baypass -gfile geno_BayPass.txt -efile env_BayPass.txt \
          -omegafile core_run_mat_omega.out \
          -outprefix aux_run -auxmodel -nthreads 4

# C2 (contrast between two pre-defined groups)
# contrasts.txt: 1 = group 1, -1 = group 2, 0 = excluded
g_baypass -gfile geno_BayPass.txt -contrastfile contrasts.txt \
          -omegafile core_run_mat_omega.out \
          -outprefix c2_run -nthreads 4

# Bayes Factor > 20 dB: strong evidence
# Bayes Factor > 30 dB: decisive evidence
# Convert: BFis (decibel scale) = 10 * log10(BF)
```

## Genomic Offset — Lind & Lotterhos 2025 Three-Regime Framework

**Goal:** Predict maladaptation under future climate using genomic-offset methods (gradientForest, RDA-offset, LFMM2offset, RONA) while characterizing prediction-novelty. For the broader genomic-prediction review (foundation of the offset literature), see Capblancq et al. 2020 *Annu Rev Ecol Evol Syst* 51:245-269.

**Approach:** Compute offset per location as the distance in transformed genomic-environmental space between current and future predicted climate. CRITICAL: characterize whether each prediction falls into regime 1 (interpolation), regime 2 (modest extrapolation), or regime 3 (highly novel, where offset is uninformative or misleading). Report offset alongside prediction-novelty mapping.

```r
library(gradientForest)
library(terra)

# Fit gradient forests on candidate loci x environment
gf <- gradientForest(cbind(env_predictors, allele_data),
                     predictor.vars = colnames(env_predictors),
                     response.vars = colnames(allele_data),
                     ntree = 500, trace = FALSE)

# Variable importance (which env drives most allele turnover)
plot(gf, plot.type = 'Overall.Importance')

# Predict genetic offset under current vs future
current <- rast('bioclim_current.tif')
future <- rast('bioclim_2070_ssp585.tif')
current_vals <- extract(current, sampling_coords)
future_vals <- extract(future, sampling_coords)

# Genomic offset = Euclidean distance in transformed space
current_transformed <- predict(gf, current_vals)
future_transformed <- predict(gf, future_vals)
genetic_offset <- sqrt(rowSums((current_transformed - future_transformed)^2))

# CRITICAL: characterize prediction-novelty regime per location
# For each future location, compute distance to nearest training-data envelope
# locations in regime 3 (high novelty) should have offset values down-weighted
# in interpretation per Lind & Lotterhos 2025

# Cross-validation with multiple offset methods (gradientForest + RDA + LFMM2offset)
# Robust signals appear in ALL methods; method-specific offset is suspect
```

## Per-Method Failure Modes

### LFMM2 with wrong K silently invalidates results

**Trigger:** Running LFMM2 with K too small (residual structure confounds environment) or K too large (latent factors absorb the environmental signal).

**Mechanism:** LFMM2 uses K latent factors to capture unobserved population structure. Wrong K means structure is either undercorrected (false positives) or overcorrected (false negatives), with no error message.

**Symptom:** Thousands of significant SNPs (too small K) OR almost none (too large K); genomic inflation factor lambda far from 1.0.

**Fix:** Use sNMF cross-entropy elbow as the primary K choice; run LFMM2 at K-1, K, K+1 and report sensitivity; flag loci detected at only one K as lower-confidence.

### RDA on un-imputed genotypes silently degrades

**Trigger:** Running RDA on an allele-frequency matrix with NA values without imputation.

**Mechanism:** RDA cannot handle missing data internally; default behavior is to use `na.omit` which drops loci AND samples. Forester 2018's RDA-dominance benchmark assumed imputed genotypes.

**Symptom:** Effective sample size much smaller than expected; results sensitive to which loci have missing data.

**Fix:** Impute missing genotypes (mean within population, kNN, or `snmf::impute()`) BEFORE running RDA.

### Genomic offset reported without prediction-novelty characterization

**Trigger:** Producing a country-wide genomic-offset map without indicating which areas are inside vs outside the training-data envelope.

**Mechanism:** Per Lind & Lotterhos 2025, genomic-offset accuracy declines with novelty of the predicted environment. Highly novel future climates fall in "regime 3" where offset is uninformative.

**Symptom:** Map shows highest predicted maladaptation in geographically extreme areas; cross-validation with multiple methods shows divergent predictions.

**Fix:** Compute prediction novelty per location (e.g., Mahalanobis distance to training envelope); annotate offset values with their regime; report cross-method consensus rather than single-method offset.

### OutFLANK applied to continuous sampling design

**Trigger:** Running OutFLANK on continuous-gradient samples treated as a single "population".

**Mechanism:** OutFLANK assumes distinct, well-defined populations for FST calculation. Continuous sampling does not have natural population delimitation; the FST distribution depends arbitrarily on how samples are binned.

**Symptom:** OutFLANK results unstable to alternative population definitions; very few or impossibly many outliers.

**Fix:** Use LFMM2 or RDA for continuous gradient data; OutFLANK is appropriate when distinct populations are defined a priori.

### LDNe physical-linkage trap on RAD-seq SNPs

**Trigger:** Running LFMM2 / RDA / OutFLANK on RAD-seq SNPs without LD pruning.

**Mechanism:** Physical linkage among RAD-seq SNPs creates correlated test statistics; outlier counts are inflated by chromosomal clustering rather than independent adaptation signal.

**Symptom:** Outliers cluster in large blocks along the genome; individual-locus inference is dominated by linkage.

**Fix:** LD-prune SNPs (PLINK `--indep-pairwise 50 5 0.2`) before GEA; OR thin to >= 1 cM apart with a genetic map.

## Quantitative Thresholds

| Threshold | Value | Source / rationale |
|-----------|-------|-------------------|
| LFMM2 K selection | Cross-entropy elbow from sNMF (`LEA`) | Caye 2019; sensitivity check at K-1, K+1 |
| Genomic inflation factor target | Lambda ~ 1.0 | > 1.5 = under-correction; < 0.5 = over-correction |
| FDR threshold | q < 0.05 (Storey) | qvalue package; more powerful than BH for genomic data |
| RDA outlier threshold | abs(z) > 3 on RDA axes 1-3 | Capblancq & Forester 2021 conservative; 2.5 for polygenic |
| BayPass Bayes Factor | BF > 20 dB strong; > 30 dB decisive | Gautier 2015 convention |
| OutFLANK trim fractions | LeftTrim = RightTrim = 0.05 | Whitlock & Lotterhos 2015 default |
| OutFLANK Hmin | 0.1 | Exclude low-heterozygosity unreliable FST |
| Forester 2018 RDA imputation | Required (no missing data) | RDA performance benchmark assumed imputed genotypes |
| Lind & Lotterhos 2025 offset regimes | Regime 1 (interpolation) / 2 (modest extrapolation) / 3 (highly novel) | Characterize per location |
| Lotterhos & Whitlock 2015 sampling design | Paired contrasts > random > transects under demographic concern | Genome-scan power |

## Common errors

| Error | Cause | Solution |
|-------|-------|----------|
| LFMM2 returns thousands of significant SNPs | K too small | Increase K via sNMF cross-entropy elbow |
| LFMM2 returns no significant SNPs | K too large; signal absorbed | Decrease K |
| Lambda (GIF) >> 1.5 | Insufficient structure correction | Increase K or add structure covariate |
| RDA effective N much smaller than total | NA values in allele matrix | Impute before RDA |
| pcadapt screeplot shows no elbow | Continuous population structure or single-cline data | Switch to LFMM2 or RDA |
| OutFLANK error about populations | Continuous sampling, no defined populations | Use LFMM2 / RDA instead |
| BayPass "format error" | Passing VCF directly | Convert to space-separated allele counts via vcf2baypass.pl |
| Gradient-forest plot blank | Insufficient SNPs in candidate set | Increase candidate-locus pool |
| Genomic offset highest in intermediate climate | PC axes capture environment; structure absorbing signal | pRDA conditioning on geographic distance |

## References

- Frichot E, Schoville SD, Bouchard G, François O (2013) LFMM. *Mol Biol Evol* 30(7):1687-1699. doi:10.1093/molbev/mst063
- Caye K, Jumentier B, Lepeule J, François O (2019) LFMM 2: fast LSE estimator. *Mol Biol Evol* 36(4):852-860. doi:10.1093/molbev/msz008
- Gain C, François O (2021) LEA 3 package. *Mol Ecol Resour* 21(8):2738-2748. doi:10.1111/1755-0998.13366
- Gautier M (2015) BayPass. *Genetics* 201(4):1555-1579. doi:10.1534/genetics.115.181453
- Forester BR, Lasky JR, Wagner HH, Urban DL (2018) Multilocus adaptation methods: RDA superior. *Mol Ecol* 27(9):2215-2233. doi:10.1111/mec.14584
- Capblancq T, Forester BR (2021) RDA Swiss-army-knife for landscape genomics. *Methods Ecol Evol* 12(12):2298-2309. doi:10.1111/2041-210X.13722
- Capblancq T, Fitzpatrick MC, Bay RA, Exposito-Alonso M, Keller SR (2020) Genomic prediction of (mal)adaptation. *Annu Rev Ecol Evol Syst* 51:245-269. doi:10.1146/annurev-ecolsys-020720-042553
- Whitlock MC, Lotterhos KE (2015) OutFLANK. *Am Nat* 186(S1):S24-S36. doi:10.1086/682949
- Lotterhos KE, Whitlock MC (2014) Demographic confounders of FST outlier tests. *Mol Ecol* 23(9):2178-2192. doi:10.1111/mec.12725
- Lotterhos KE, Whitlock MC (2015) Sampling-design optimization for genome scans. *Mol Ecol* 24(5):1031-1046. doi:10.1111/mec.13100
- Lind BM, Lotterhos KE (2025) Accuracy of predicting maladaptation. *Mol Ecol Resour* 25:e14008. doi:10.1111/1755-0998.14008
- Wang IJ, Bradburd GS (2014) Isolation by environment (IBD vs IBE). *Mol Ecol* 23(23):5649-5662. doi:10.1111/mec.12938
- McRae BH, Dickson BG, Keitt TH, Shah VB (2008) Circuitscape circuit-theory. *Ecology* 89(10):2712-2724. doi:10.1890/07-1861.1
- Peterman WE (2018) ResistanceGA. *Methods Ecol Evol* 9(6):1638-1647. doi:10.1111/2041-210X.12984
- Ellis N, Smith SJ, Pitcher CR (2012) Gradient forests. *Ecology* 93(1):156-168. doi:10.1890/11-0252.1
- Aguirre-Liguori JA, Ramírez-Barahona S, Gaut BS (2021) Evolutionary genomics of climate response. *Nat Ecol Evol* 5(10):1350-1360. doi:10.1038/s41559-021-01526-9
- Hoban S, Kelley JL, Lotterhos KE et al. (2016) Finding the genomic basis of local adaptation. *Am Nat* 188(4):379-397. doi:10.1086/688018
- Luu K, Bazin E, Blum MGB (2017) pcadapt: an R package to perform genome scans for selection based on principal component analysis. *Mol Ecol Resour* 17(1):67-77. doi:10.1111/1755-0998.12592
- Legendre P, Fortin M-J (2010) Comparison of the Mantel test and alternative approaches for detecting complex multivariate relationships. *Mol Ecol Resour* 10(5):831-844. doi:10.1111/j.1755-0998.2010.02866.x

## Related Skills

- ecological-genomics/conservation-genetics - Population genetic health (Ne, F_ROH); use ESU/MU framework not species delimitation for management units
- ecological-genomics/community-ecology - Environmental gradient analysis for species composition (analog of GEA for community data)
- population-genetics/selection-statistics - Selection scans in human population genetics
- population-genetics/population-structure - STRUCTURE / ADMIXTURE for substructure inference
- variant-calling/vcf-basics - VCF preparation from RAD-seq or WGS
