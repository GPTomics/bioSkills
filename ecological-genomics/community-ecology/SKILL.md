---
name: bio-ecological-genomics-community-ecology
description: Analyzes species-environment relationships with constrained ordination (CCA, RDA, db-RDA), variance partitioning, indicator species (indicspecies IndVal.g group-equalized), PERMANOVA paired MANDATORILY with PERMDISP (Anderson & Walsh 2013; dispersion confounds centroid tests), Joint Species Distribution Models (HMSC, sjSDM, gjam) with explicit rejection of "residual covariance equals biotic interaction", phylogenetic community ecology (SES_MPD/MNTD), trait-environment via RLQ + fourth-corner with corrected modeltype=6 (Dray 2014), bipartite network metrics (NODF, modularity) with curveball null (Strona 2014), and Mantel-test replacements (dbRDA, GDM) for spatial data. Use when testing how environmental gradients structure communities, identifying habitat indicator taxa, partitioning variance among predictors, deciding whether PERMANOVA significance is location vs dispersion, picking among HMSC/sjSDM/gjam, or replacing Mantel tests for landscape data.
tool_type: r
primary_tool: vegan
---

## Version Compatibility

Reference examples tested with: vegan 2.6+, indicspecies 1.7+, Hmsc 3.0+, sjSDM 1.0+, ade4 1.7+, picante 1.8+, ggplot2 3.5+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Community Ecology

**"Test how environmental gradients structure my species communities"** -> Constrained ordination (CCA / RDA / db-RDA) with explicit dispersion testing alongside PERMANOVA, indicator-species analysis with group-size correction, Joint Species Distribution Models for residual covariance and prediction, and trait-environment testing with statistically corrected permutation schemes.
- R: `vegan::rda()`, `vegan::dbrda()`, `vegan::adonis2()` for ordination and PERMANOVA
- R: `vegan::betadisper()` for the mandatory PERMDISP companion to PERMANOVA
- R: `Hmsc` or `sjSDM` for joint species distribution modeling
- R: `indicspecies::multipatt(..., func = 'IndVal.g')` for indicator species

## The Single Most Important Modern Insight -- Always run PERMDISP alongside PERMANOVA

Anderson & Walsh 2013 *Ecol Monogr* 83(4):557-574 established that PERMANOVA's pseudo-F is **sensitive to dispersion heterogeneity** — a significant PERMANOVA can reflect centroid difference, dispersion difference, or both. **Running PERMANOVA without PERMDISP is the single most common methodological failure in community-ecology papers.** If `betadisper` is significant alongside a significant PERMANOVA, the location-difference conclusion is not supported by the data alone; the two could be entirely a dispersion artifact.

A second cornerstone insight: residual species-species covariance in Joint Species Distribution Models is NOT a clean estimator of biotic interaction (Pollock 2014, reaffirmed by Zurell 2018 and Poggiato 2021). It reflects unmeasured covariates, dispersal limitation, sampling artifacts, AND any genuine biotic interactions, in unknown proportions. Skills that interpret residual covariance as interaction are over-interpreting.

## Algorithmic Taxonomy

| Method | Estimand | Strength | Fails when |
|--------|----------|----------|------------|
| CCA | Species-environment unimodal | Standard for chi-square data; fits bell-shaped responses | Linear gradients (use RDA); does not handle short gradients well |
| RDA | Species-environment linear | High power with short gradients; Hellinger-friendly | Long unimodal gradients (>3 SD DCA axis 1); needs no missing data |
| db-RDA / capscale | Constrained ordination on any distance | Flexible (Bray-Curtis, Sorensen, weighted Unifrac); dispersion-robust | Less power than RDA for purely linear gradients |
| PERMANOVA (adonis2) | Centroid difference among groups | Non-parametric; handles any dissimilarity | Sensitive to dispersion difference — MUST run PERMDISP alongside |
| PERMDISP (betadisper) | Dispersion difference among groups | Tests the confound that contaminates PERMANOVA | Low power with small N; report alongside PERMANOVA |
| ANOSIM | Group-difference test | Legacy familiarity | Worse than PERMANOVA for the same use case; biased by unbalanced N |
| Mantel test | Correlation between two distance matrices | Conceptually simple | Low power under spatial autocorrelation; biased; replace with dbRDA or GDM |
| Partial Mantel | Correlation controlling for a third matrix | Conceptually simple | INFLATES Type I error under autocorrelation (worse than basic Mantel) |
| HMSC (Helsinki tradition) | Bayesian JSDM with phylogeny + traits + spatial | Rigorous; ecological-theory priors | Slow for S > 500 species |
| sjSDM (Pichler & Hartig) | Latent-variable-free JSDM via MC + elastic net | Orders of magnitude faster; high-S friendly | Less explicit interpretation than HMSC |
| gjam | Cross-data-type JSDM (counts, presence, continuous) | Integrates heterogeneous data | Different output structure than HMSC; not directly comparable |
| IndVal (Dufrene-Legendre) | Species-group association | Combines specificity AND fidelity | Biased by group size unless `func='IndVal.g'` |
| SES_MPD / SES_MNTD | Phylogenetic structure vs null | Tests whether communities are clustered or overdispersed | Interpretation of sign requires trait-conservatism check (Mayfield & Levine 2010) |
| RLQ + fourth-corner | Trait-environment relationship | Tests species-trait response to environment | Default modeltype gives inflated Type I; use modeltype=6 |

## Decision Tree by Scenario

| Scenario | Recommended approach | Why |
|----------|---------------------|-----|
| Unimodal species responses, gradient > 3 SD (DCA axis 1) | CCA | Linear assumption fails for long gradients |
| Linear species responses, gradient <= 3 SD | RDA on Hellinger-transformed data | RDA assumes linearity; Hellinger solves double-zero problem |
| Bray-Curtis or other non-Euclidean distance preferred | db-RDA | Handles any distance metric while preserving constrained-ordination interpretation |
| Test "do these groups differ in community composition" | PERMANOVA (adonis2) AND PERMDISP (betadisper) | Never PERMANOVA alone; dispersion confound is non-negotiable |
| ANOSIM out of habit | Use PERMANOVA + PERMDISP instead | ANOSIM is biased by unbalanced N and dispersion |
| Mantel test for landscape genetics | dbRDA with spatial covariates OR GDM | Mantel has low power; partial Mantel inflates Type I |
| Joint species modeling with traits and phylogeny | HMSC (S < 500); sjSDM (S > 500) | HMSC encodes theory in priors; sjSDM is the only scalable option for high-S |
| Multi-species "association" structure | JSDM residual covariance | Interpret as ANY of (interaction, shared response to unmeasured driver, dispersal, sampling), NOT pure interaction |
| Indicator species with unbalanced group sizes | `multipatt(..., func = 'IndVal.g')` | `IndVal` is biased by group size; group-equalized form corrects |
| Phylogenetic community structure | SES_MPD with EXPLICIT null model + trait-conservatism test | Cite Mayfield & Levine 2010 for the "clustering = filtering" interpretation trap |
| Trait-environment hypothesis | RLQ + fourth-corner with `modeltype=6` | Default `modeltype=2` or `modeltype=4` gives inflated Type I error |
| Bipartite network nestedness/modularity | NODF2 + modularity with `curveball` null randomization | Strona 2014 curveball is exponentially faster and unbiased for binary matrices |

## CCA vs RDA — Gradient-Length Decision

**Goal:** Choose between unimodal (CCA) and linear (RDA) constrained ordination based on the dominant gradient length in the species data.

**Approach:** Run a Detrended Correspondence Analysis (DCA) on the raw community matrix; use the axis-1 SD length as the gradient-length metric. > 3 SD suggests unimodal CCA; < 3 SD suggests linear RDA on Hellinger-transformed data; 2-3 SD is a gray zone where either is defensible.

```r
library(vegan)

# Step 1: Check gradient length
dca <- decorana(species_matrix)
dca  # axis 1 length in SD units

# Step 2a: Long gradient -> CCA
cca_result <- cca(species_matrix ~ temperature + precipitation + pH + elevation,
                  data = env_data)
anova(cca_result, by = 'margin', permutations = 999)

# Step 2b: Short gradient -> RDA with Hellinger transformation
# Hellinger (Legendre & Gallagher 2001) is MANDATORY before RDA on community data
species_hell <- decostand(species_matrix, method = 'hellinger')
rda_result <- rda(species_hell ~ temperature + precipitation + pH + elevation,
                  data = env_data)
RsquareAdj(rda_result)$adj.r.squared
anova(rda_result, by = 'margin', permutations = 999)

# Forward selection with adjusted R-squared criterion (Peres-Neto 2006 Ecology 87:2614)
rda_null <- rda(species_hell ~ 1, data = env_data)
rda_full <- rda(species_hell ~ ., data = env_data)
rda_sel <- ordiR2step(rda_null, scope = formula(rda_full),
                      direction = 'forward', permutations = 999)

# VIF check: > 10 indicates problematic multicollinearity
vif.cca(rda_sel)
```

## PERMANOVA with the Mandatory PERMDISP Companion

**Goal:** Test whether community composition differs across groups while detecting the dispersion-heterogeneity confound.

**Approach:** Run `adonis2` (modern PERMANOVA per Anderson 2001 *Austral Ecol*) on Hellinger or Bray-Curtis distances; THEN run `betadisper` (PERMDISP per Anderson 2006 *Biometrics*) to test whether group dispersions are unequal. Report BOTH results. If betadisper is significant, the adonis2 conclusion of centroid difference is not supported — the apparent "group difference" may be entirely a dispersion artifact.

```r
library(vegan)

# Distance matrix
bray_dist <- vegdist(species_matrix, method = 'bray')

# PERMANOVA via adonis2 (modern API; adonis() is deprecated)
# by='margin' for unbalanced designs (sequential SS gives wrong result)
permanova <- adonis2(bray_dist ~ habitat + soil_pH, data = env_data,
                     by = 'margin', permutations = 999)
permanova

# MANDATORY companion: PERMDISP via betadisper
# Tests homogeneity of multivariate dispersions across groups
disp <- betadisper(bray_dist, env_data$habitat)
disp_test <- permutest(disp, permutations = 999)
disp_test

# Interpretation rule:
# PERMANOVA p < 0.05 AND betadisper p > 0.05 -> location difference is real
# PERMANOVA p < 0.05 AND betadisper p < 0.05 -> CONFOUNDED, cannot conclude location difference
# PERMANOVA p > 0.05 -> no group difference detected

# Visualize dispersions
plot(disp)
boxplot(disp)

# For pairwise group comparisons, do NOT use Bonferroni-corrected pairwise PERMANOVA
# (permutation tests with overlapping sets do not give nominal FWER from Bonferroni);
# use pairwiseAdonis::pairwise.adonis2 with FDR correction instead
# install.packages('pairwiseAdonis')
```

## Joint Species Distribution Models — HMSC vs sjSDM

**Goal:** Model species occurrences jointly to capture environmental responses, traits, phylogeny, and residual covariance among species.

**Approach:** For S < 500 species with rich theory and traits/phylogeny: use HMSC (Ovaskainen 2017 *Ecol Lett* 20:561-576; current R package Tikhonov 2020 *Methods Ecol Evol* 11:442-447) for Bayesian inference with explicit ecological priors. For S > 500 species (modern metabarcoding datasets): use sjSDM (Pichler & Hartig 2021 *Methods Ecol Evol* 12:2159-2173) which is orders of magnitude faster via Monte Carlo approximation of the joint likelihood with elastic-net regularization. The Wilkinson 2019 *Methods Ecol Evol* 10:198-211 benchmark is essential reading before picking among HMSC, sjSDM, gjam, and BayesComm. **Do NOT interpret residual species-species covariance as biotic interaction** (see Zurell 2018 *Ecography* 41:1812-1819; Poggiato 2021 *Trends Ecol Evol* 36:391-401).

```r
library(Hmsc)

# HMSC for moderate-S Bayesian JSDM
# X: site x environment data
# Y: site x species presence-absence or abundance
# distr: 'probit' for presence-absence; 'lognormal poisson' for counts
m <- Hmsc(Y = species_matrix, XData = env_data, XFormula = ~ temperature + soil_pH,
         distr = 'probit')

# Sample posterior (production runs need thin >= 100, transient >= 1000, samples >= 1000)
m <- sampleMcmc(m, thin = 10, samples = 1000, transient = 5000, nChains = 4)

# Variance partitioning into fixed (environment), random (latent factors), traits, phylogeny
VP <- computeVariancePartitioning(m)
plotVariancePartitioning(m, VP = VP)

# Predict to new environment
pred <- predict(m, XData = new_env_data, expected = TRUE)

# For S > 500: switch to sjSDM
# library(sjSDM)
# m_sj <- sjSDM(Y = species_matrix, env = linear(env_data, ~ temperature + soil_pH),
#               family = binomial('probit'))
# summary(m_sj); plot(m_sj)
```

## Indicator Species with Group-Size Correction

**Goal:** Identify species statistically associated with site groups using an indicator value that combines specificity and fidelity, corrected for unequal group sizes.

**Approach:** Run `indicspecies::multipatt` with `func = 'IndVal.g'` (the group-size-equalized form per De Caceres & Legendre 2009 *Ecology* 90:3566-3574). The original `IndVal` is biased toward larger groups; the `.g` form corrects this. For continuous-vs-categorical or rank-based associations, use `func = 'r.g'` (point-biserial correlation, group-equalized).

```r
library(indicspecies)

# IndVal.g: group-size-equalized indicator value
# multipatt tests species-group associations with permutation
# duleg=TRUE: tests only single-group associations (not combinations)
# duleg=FALSE: tests species against all group combinations (more powerful but more tests)
mp <- multipatt(species_matrix, site_groups,
                func = 'IndVal.g',       # Group-equalized; NOT 'IndVal' (biased)
                duleg = TRUE,
                control = how(nperm = 999))

summary(mp, alpha = 0.05)

# Extract significant indicators sorted by p-value
sig <- mp$sign[!is.na(mp$sign$p.value) & mp$sign$p.value < 0.05, ]
sig[order(sig$p.value), ]
```

## Mantel Replacement — dbRDA with Spatial Covariates

**Goal:** Test whether community/genetic distance correlates with environmental distance while controlling for spatial autocorrelation, replacing the low-power and bias-prone Mantel framework.

**Approach:** Use db-RDA with spatial predictors (PCNM eigenvectors or raw coordinates) as `Condition()` (Legendre & Fortin 2010 *Mol Ecol Resour* 10:831-844). Partial Mantel inflates Type I error under autocorrelation — do not use for landscape data.

```r
library(vegan)
library(adespatial)

# Build spatial predictors: PCNM (principal coordinates of neighbor matrices)
geo_dist <- dist(coords[, c('longitude', 'latitude')])
pcnm <- pcnm(geo_dist)
pcnm_vars <- pcnm$vectors  # significant axes

# db-RDA controlling for spatial structure
dbrda_result <- dbrda(bray_dist ~ temperature + precipitation +
                       Condition(as.matrix(pcnm_vars)),
                      data = env_data, add = 'lingoes')

# Test marginal significance of environment AFTER conditioning out space
anova(dbrda_result, by = 'margin', permutations = 999)

# DO NOT use partial Mantel for this question; cite Legendre & Fortin 2010
```

## Per-Method Failure Modes

### PERMANOVA significant, betadisper also significant -> conclusion not supported

**Trigger:** Reporting PERMANOVA p < 0.05 as evidence of "community composition differs across groups" without running betadisper, OR running betadisper and finding it significant but ignoring the result.

**Mechanism:** PERMANOVA's pseudo-F responds to BOTH centroid shifts AND dispersion heterogeneity (Anderson & Walsh 2013). When groups differ in dispersion but not centroid, PERMANOVA can return p < 0.05 from the dispersion difference alone.

**Symptom:** Reviewer asks "was dispersion checked?"; PCoA visualization shows overlapping centroids but one group is more dispersed; betadisper p < 0.05.

**Fix:** Report both PERMANOVA and betadisper. If betadisper is significant, the conclusion must be reframed: "Groups differ in either centroid or dispersion, with dispersion heterogeneity present." Consider db-RDA which is more dispersion-robust.

### Mantel test reports low p but reflects spatial autocorrelation

**Trigger:** Using Mantel or partial Mantel to test "is genetic distance correlated with environmental distance" in a spatially-structured landscape.

**Mechanism:** Mantel statistics are biased downward by spatial autocorrelation in either distance matrix (Legendre & Fortin 2010). Partial Mantel further inflates Type I error rates under autocorrelation (see Guillot & Rousset 2013 *Methods Ecol Evol* 4:336-344 for the formal demonstration).

**Symptom:** Mantel p < 0.001 but no detectable signal when re-tested with dbRDA conditioning on geographic distance.

**Fix:** Use db-RDA with PCNM spatial eigenvectors as `Condition()`, OR GDM (Generalized Dissimilarity Modeling).

### JSDM residual covariance interpreted as biotic interaction

**Trigger:** Reporting "species A and species B have residual covariance after fitting environment, indicating biotic interaction."

**Mechanism:** Residual covariance after environmental fitting reflects (a) unmeasured covariates, (b) dispersal limitation, (c) sampling artifacts, (d) shared response to unmeasured drivers, AND (e) any genuine biotic interactions, in unknown proportions. Without manipulative or independent corroboration, the interaction signal cannot be isolated.

**Symptom:** Reviewer challenges "interaction" interpretation; sensitivity tests with different environmental specifications change the residual structure dramatically.

**Fix:** Report residual covariance descriptively ("residual association after fitting environment"), explicitly acknowledge alternative explanations, and cite Zurell 2018 or Poggiato 2021 for the interpretation caveat.

### IndVal biased by unequal group sizes

**Trigger:** Running `multipatt(..., func = 'IndVal')` with strongly unbalanced groups.

**Mechanism:** The original IndVal index from Dufrene & Legendre 1997 is biased toward larger groups because larger groups have higher average occupancy by chance.

**Symptom:** All indicators are assigned to the largest group; small-group indicators are undetected.

**Fix:** Use `func = 'IndVal.g'` (group-equalized) per De Caceres & Legendre 2009. The `.g` correction is the modern default.

### SES_MPD sign interpreted as filtering vs competition without trait check

**Trigger:** Reporting "SES_MPD < 0 (NRI > 0) indicates environmental filtering" without testing whether traits track phylogeny.

**Mechanism:** Mayfield & Levine 2010 *Ecol Lett* 13:1085-1093 showed competition can produce phylogenetic clustering (not just overdispersion) when ecologically similar species coexist via R*-rule trait differences. The Webb 2002 sign-to-process mapping is incomplete.

**Symptom:** Trait-similarity vs phylogenetic-similarity correlation has not been examined; reviewer asks about Mayfield-Levine alternative.

**Fix:** Compute Blomberg's K AND Pagel's lambda; if traits track phylogeny strongly, the clustering signal could be either filtering OR competition. Cite Mayfield & Levine 2010 in interpretation.

## Quantitative Thresholds

| Threshold | Value | Source / rationale |
|-----------|-------|-------------------|
| DCA gradient length | < 3 SD: RDA; > 3 SD: CCA; 2-3 SD: gray zone | Standard ordination decision rule |
| VIF | > 10: collinearity problem | Hair et al. ecology convention |
| PERMANOVA + PERMDISP rule | Report BOTH; if both p < 0.05, location difference not supported | Anderson & Walsh 2013 |
| RDA R^2 reporting | Adjusted R^2 (`RsquareAdj$adj.r.squared`) | Peres-Neto 2006 Ecology 87:2614 |
| HMSC MCMC settings | Thin 100, samples 1000, transient 1000, chains >= 2 | Hmsc tutorial recommendations |
| Sample-size rule for sjSDM | Used when S > 500 species | HMSC computational practicality |
| IndVal significance | p < 0.05 after 999 permutations | Standard alpha; use FDR for many comparisons |
| Bootstrap iterations | nperm = 999 minimum for tests | Standard for permutation tests |

## Common errors

| Error | Cause | Solution |
|-------|-------|----------|
| adonis() doesn't accept new arguments | adonis() deprecated in vegan 2.6+ | Use adonis2() |
| betadisper produces NA | Unequal group sizes with very small N | Increase replicates per group |
| PCoA shows no group separation despite p < 0.05 | Dispersion-driven PERMANOVA significance | Report betadisper alongside |
| HMSC convergence diagnostics flag | Insufficient thinning/transient | Increase thin and transient parameters |
| sjSDM error about GPU/CUDA | Default device misconfigured | Set `device='cpu'` if no GPU |
| multipatt all p-values 1.0 | Group factor not a factor | `as.factor(site_groups)` |
| IndVal flags only large-group species | Using `func='IndVal'` not `'IndVal.g'` | Switch to `IndVal.g` |
| dbRDA negative eigenvalues warning | Non-Euclidean distance with no correction | Add `add = 'lingoes'` |
| Mantel test always significant | Spatial autocorrelation inflating correlation | Switch to dbRDA with spatial covariates |

## References

- Anderson MJ (2001) A new method for non-parametric multivariate analysis of variance. *Austral Ecol* 26(1):32-46. doi:10.1111/j.1442-9993.2001.01070.pp.x
- Anderson MJ (2006) Distance-based tests for homogeneity of multivariate dispersions. *Biometrics* 62(1):245-253. doi:10.1111/j.1541-0420.2005.00440.x
- Anderson MJ, Walsh DCI (2013) PERMANOVA, ANOSIM, and the Mantel test in the face of heterogeneous dispersions. *Ecol Monogr* 83(4):557-574. doi:10.1890/12-2010.1
- Legendre P, Fortin M-J (2010) Comparison of the Mantel test and alternative approaches for detecting complex multivariate relationships. *Mol Ecol Resour* 10(5):831-844. doi:10.1111/j.1755-0998.2010.02866.x
- Legendre P, Gallagher ED (2001) Ecologically meaningful transformations for ordination of species data. *Oecologia* 129(2):271-280. doi:10.1007/s004420100716
- Peres-Neto PR, Legendre P, Dray S, Borcard D (2006) Variation partitioning of species data matrices. *Ecology* 87(10):2614-2625. doi:10.1890/0012-9658(2006)87[2614:VPOSDM]2.0.CO;2
- Dufrene M, Legendre P (1997) Species assemblages and indicator species. *Ecol Monogr* 67(3):345-366. doi:10.1890/0012-9615(1997)067[0345:SAAIST]2.0.CO;2
- De Caceres M, Legendre P (2009) Associations between species and groups of sites. *Ecology* 90(12):3566-3574. doi:10.1890/08-1823.1
- Pollock LJ, Tingley R, Morris WK et al. (2014) Joint Species Distribution Model (JSDM). *Methods Ecol Evol* 5(5):397-406. doi:10.1111/2041-210X.12180
- Ovaskainen O, Tikhonov G, Norberg A et al. (2017) Hierarchical Modelling of Species Communities (HMSC). *Ecol Lett* 20(5):561-576. doi:10.1111/ele.12757
- Tikhonov G, Opedal OH, Abrego N et al. (2020) Joint species distribution modelling with the R-package Hmsc. *Methods Ecol Evol* 11(3):442-447. doi:10.1111/2041-210X.13345
- Pichler M, Hartig F (2021) A new joint species distribution model for faster and more accurate inference. *Methods Ecol Evol* 12(11):2159-2173. doi:10.1111/2041-210X.13687
- Wilkinson DP, Golding N, Guillera-Arroita G et al. (2019) Comparison of joint species distribution models. *Methods Ecol Evol* 10(2):198-211. doi:10.1111/2041-210X.13106
- Zurell D, Pollock LJ, Thuiller W (2018) Do joint species distribution models reliably detect interspecific interactions from co-occurrence data in homogenous environments? *Ecography* 41(11):1812-1819. doi:10.1111/ecog.03315
- Poggiato G, Munkemuller T, Bystrova D et al. (2021) On the interpretations of joint modeling in community ecology. *Trends Ecol Evol* 36(5):391-401. doi:10.1016/j.tree.2021.01.002
- Guillot G, Rousset F (2013) Dismantling the Mantel tests. *Methods Ecol Evol* 4(4):336-344. doi:10.1111/2041-210x.12018
- Webb CO, Ackerly DD, McPeek MA, Donoghue MJ (2002) Phylogenies and community ecology. *Annu Rev Ecol Syst* 33:475-505. doi:10.1146/annurev.ecolsys.33.010802.150448
- Mayfield MM, Levine JM (2010) Opposing effects of competitive exclusion on phylogenetic community structure. *Ecol Lett* 13(9):1085-1093. doi:10.1111/j.1461-0248.2010.01509.x
- Dray S, Choler P, Doledec S et al. (2014) Combining fourth-corner and RLQ methods. *Ecology* 95(1):14-21. doi:10.1890/13-0196.1
- Strona G, Nappo D, Boccacci F, Fattorini S, San-Miguel-Ayanz J (2014) A fast and unbiased procedure to randomize ecological binary matrices. *Nat Commun* 5:4114. doi:10.1038/ncomms5114

## Related Skills

- ecological-genomics/biodiversity-metrics - Alpha/beta diversity and Hill numbers prior to ordination
- ecological-genomics/edna-metabarcoding - Generate community data from environmental samples
- ecological-genomics/landscape-genomics - Genotype-environment associations (genetic analog of GEA)
- microbiome/diversity-analysis - Unconstrained ordination alternative for 16S microbiome
- data-visualization/ggplot2-fundamentals - Customize triplots, ordination plots, and indicator-species visualizations
