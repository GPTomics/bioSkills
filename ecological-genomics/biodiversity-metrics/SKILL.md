---
name: bio-ecological-genomics-biodiversity-metrics
description: Quantifies biodiversity from species abundance/incidence tables using Hill numbers (iNEXT) with coverage-based rarefaction-extrapolation (Chao & Jost 2012), asymptotic richness via Chao1/ACE/jackknife as a lower bound, Baselga turnover/nestedness partition with the Podani alternative as sensitivity check, mandatory Hellinger transformation before ordination (Legendre & Gallagher 2001), Faith PD and SES_MPD/SES_MNTD with explicit null-model choice, and Maire 2015 functional-diversity dimensionality optimization. Use when comparing diversity across sites with unequal sampling effort, picking the right richness estimator for singleton-heavy amplicon data, partitioning beta diversity into turnover vs nestedness, reporting Hill-number effective species counts rather than raw entropies, computing SES_MPD with explicit null-model justification, or deciding whether to apply standard metrics to compositional amplicon data. Not for clinical 16S microbiome diversity (see microbiome/diversity-analysis).
tool_type: r
primary_tool: iNEXT
---

## Version Compatibility

Reference examples tested with: iNEXT 3.0+, iNEXT.3D 1.0+, vegan 2.6+, betapart 1.6+, picante 1.8+, ggplot2 3.5+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Biodiversity Metrics

**"Calculate species diversity for my ecological samples"** -> Compute Hill-number diversity (numbers-equivalent of richness, Shannon, Simpson) with coverage-based rarefaction/extrapolation, choose the richness estimator appropriate to the singleton/doubleton signature of the data, and partition beta diversity into turnover and nestedness components with a documented partition framework.
- R: `iNEXT::iNEXT()` for coverage-based rarefaction/extrapolation
- R: `betapart::beta.multi()` for Baselga turnover/nestedness partition
- R: `picante::ses.mpd()` for phylogenetic-community SES with explicit null model

## The Single Most Important Modern Insight -- Standardize by COVERAGE not by sample size

Chao & Jost 2012 *Ecology* 93(12):2533-2547 established that comparing diversity across assemblages by rarefying to a common sample size systematically biases comparisons whenever assemblages differ in underlying diversity: a 100-read rarefaction of a 50-species community is essentially saturated (coverage approximately 99%) while the same 100 reads from a 500-species community covers only approximately 40% of the underlying diversity. The two rarefied diversities are not measuring the same thing. **Coverage-based rarefaction with iNEXT is the postdoc-grade default; sample-size rarefaction is now considered a methodological anti-pattern for cross-site comparison.**

A second insight pairs with this: raw Shannon and Simpson indices are entropies, NOT diversities. Jost 2006 *Oikos* 113(2):363-375 showed that only their numbers-equivalents (exp(H), 1/D) are comparable as effective species counts and have the intuitive "doubling property" (merging two equally diverse equally abundant assemblages doubles the diversity). Reporting raw Shannon = 3.2 vs 3.0 hides whether the difference is large or trivial; reporting `1`D = 24.5 vs 20.1 species-equivalents makes the 22% gap visible.

## Algorithmic Taxonomy

| Method | Estimand | Strength | Fails when |
|--------|----------|----------|------------|
| Hill numbers q=0,1,2 (iNEXT) | Effective species count at order q | Unifies richness, Shannon, Simpson; doubling property | None; report all three q values |
| Chao1 = S + f1^2/(2*f2) | Lower bound on asymptotic richness | Non-parametric; works at low coverage | Singletons dominated by PCR error (amplicon data); f2 near zero |
| Chao1bc = S + f1(f1-1)/(2(f2+1)) | Bias-corrected Chao1 | Stable when f2 = 0 | Same singleton-bias issue |
| ACE | Asymptotic richness using all rare classes (f1...f10) | Uses more rare-class information than Chao1 | Choice of "rare" cutoff (default 10) is arbitrary |
| Jackknife1 = S + f1*(n-1)/n | Asymptotic richness via resampling | Robust when f2 = 0 | Sensitive to singleton count alone |
| Coverage-based rarefaction (iNEXT) | Diversity at standardized completeness | Correct comparison across sites with unequal effort | Extrapolation beyond 2x reference size is unreliable |
| Sample-size rarefaction | Diversity at common n | Legacy familiar | Systematically biased when communities differ in true diversity |
| Faith's PD | Sum of branch lengths on spanning tree | Captures evolutionary distinctness | Sensitive to richness; report alongside SES_PD |
| Rao's Q | Pairwise functional/phylogenetic distance * abundance | Unifies taxonomic and functional/phylogenetic diversity | Trait dimensionality artifacts (see Maire 2015) |
| FRic / FEve / FDiv | Functional richness/evenness/divergence | Multidimensional trait coverage | FRic inflates with collinear traits; optimize axis count |
| Baselga beta partition | Sorensen = turnover (Simpson) + nestedness | Decomposes beta into ecological processes | Not unique; Podani partition gives different components |
| Podani / Carvalho partition | Alternative turnover + richness-difference | Conceptually distinct from Baselga | Same data, different conclusion possible; report both |
| Hellinger transform + PCA/RDA | Solves double-zero problem | Standard for community ordination | None; required preprocessing |

## Decision Tree by Scenario

| Scenario | Recommended approach | Why |
|----------|---------------------|-----|
| Cross-site diversity comparison with unequal sampling effort | Coverage-based rarefaction in iNEXT to common C (typically 0.95) | Sample-size rarefaction is biased when assemblages differ in true diversity |
| Reporting diversity for publication | Hill numbers q=0,1,2 (numbers-equivalents) | Comparable, additive, doubling property; raw entropies are not |
| Amplicon/eDNA data with many singletons | Skip Chao1 OR use after careful denoising; report Good's coverage alongside | Singletons in amplicon data are dominated by PCR error, not undersampling |
| Singleton-heavy real data (f1 >> f2) | Chao1 with wide CI; cross-check with ACE | Chao1 variance grows as f1^4; ACE uses more rare-class information |
| f2 = 0 (no doubletons) | Chao1bc or jackknife1 | Original Chao1 is undefined; bias-corrected form is the standard |
| Diversity at much larger sample size than observed | Stop extrapolation at 2x reference size (the doubling rule) | iNEXT will silently extrapolate further; variance grows superlinearly beyond |
| Beta diversity for two assemblages | Sorensen or Jaccard with Baselga partition; document choice | Bray-Curtis is not a true metric and breaks some downstream methods |
| Beta diversity reported as turnover/nestedness | Run BOTH Baselga AND Podani partitions and report both | The partition is not mathematically unique; one alone hides ambiguity |
| Before PCA/RDA on community data | Hellinger transformation first | Solves double-zero problem; raw Euclidean PCA on counts is malpractice |
| Phylogenetic community structure (NRI/NTI) | SES_MPD with explicit null.model justification | Default null may not match question; document choice |
| Functional diversity (FRic) | Optimize PCoA axis count via Maire 2015 mAD; report chosen k | FRic biased by trait-axis count; defaults (k=2-3) are too few |
| Amplicon/compositional data | Either presence/absence metrics (Chao2, Sorensen) OR explicit CLR-based diversity per Gloor 2017 *Front Microbiol* 8:2224 | Hill-numbers assume absolute abundances; amplicon counts are compositional |

## Hill Numbers and Why Raw Entropies Are Not Diversities

**Goal:** Report diversity values that are interpretable as effective species counts and additive under standard partitions.

**Approach:** Compute Hill numbers `q`D at q = 0, 1, 2 (Jost 2006 *Oikos* 113:363-375; iNEXT software paper Hsieh, Ma, Chao 2016 *Methods Ecol Evol* 7:1451-1456). q controls sensitivity to rare vs common species in a continuous parametric family: q=0 weights all species equally (richness); q=1 is the geometric-mean weighting (Shannon-equivalent exp(H)); q=2 weights down rare species rapidly (Simpson-equivalent 1/D). qD has the doubling property — merging two equally diverse equally abundant assemblages exactly doubles qD. Raw Shannon and Simpson values do NOT.

```r
library(iNEXT)
library(vegan)

abundance_data <- list(
    site_A = c(100, 45, 23, 12, 8, 5, 3, 2, 1, 1),
    site_B = c(80, 60, 40, 30, 20, 15, 10, 8, 5, 3, 2, 1),
    site_C = c(200, 10, 5, 2, 1, 1, 1)
)

# Hill numbers q=0,1,2 with coverage-based rarefaction/extrapolation
# nboot=200 is the publication-quality floor; default nboot=50 is too few for CIs
result <- iNEXT(abundance_data, q = c(0, 1, 2), datatype = 'abundance', nboot = 200)

# Report effective species counts at standardized coverage (postdoc-grade default)
# coverage=0.95: 95% of individuals in the underlying community belong to detected species
est <- estimateD(abundance_data, datatype = 'abundance', base = 'coverage', level = 0.95)
est

# For vegan equivalence: numbers-equivalent of Shannon
shannon_eff <- exp(diversity(community_matrix, index = 'shannon'))
# Inverse Simpson is already in numbers-equivalent units
invsimp <- diversity(community_matrix, index = 'invsimpson')
```

## Asymptotic Richness — Chao1 is a Lower Bound, NOT a Point Estimate

**Goal:** Estimate the lower bound on true richness when sampling is incomplete and report it correctly.

**Approach:** Compute Chao1 from the singleton/doubleton ratio (Chao 1984), check that singletons are real biology rather than PCR/sequencing error, and report as a LOWER BOUND with Good's coverage to indicate sampling adequacy. Switch to ACE or jackknife1 when f2 is near zero or singletons are unreliable.

The original Chao1 derivation (Chao 1984) under a Gamma-Poisson mixture model produces a NON-PARAMETRIC LOWER BOUND on richness, not a point estimate. The bound is tight only under specific homogeneity assumptions. Reporting "Chao1 estimated richness = 320" without "at least" is widespread in published literature but statistically wrong.

```r
library(iNEXT)

# AsyEst returns Chao1 (q=0), Chao-Shannon (q=1), Chao-Simpson (q=2) with CIs
asymp <- iNEXT(abundance_data, q = c(0, 1, 2), datatype = 'abundance')$AsyEst

# CRITICAL: also report Good's coverage to indicate whether the bound is informative
# Coverage < 0.85 means heavily under-sampled; Chao1 CI will be huge and bound is uninformative
coverage <- estimateD(abundance_data, datatype = 'abundance',
                      base = 'coverage', level = 0.95)
# Inspect estimateD output's Coverage column at the OBSERVED sample size

# Singleton check: if f1 >> f2 and singletons are likely PCR artifacts (amplicon data),
# do NOT report Chao1 — it will measure "how much PCR error" not "how much undiscovered diversity"
f1_check <- sapply(abundance_data, function(x) sum(x == 1))
f2_check <- sapply(abundance_data, function(x) sum(x == 2))
cat('Singletons f1:', f1_check, '\n')
cat('Doubletons f2:', f2_check, '\n')
cat('f1^2/(2*f2) Chao1 bound term:', f1_check^2 / (2 * pmax(f2_check, 0.5)), '\n')
```

## Coverage-Based Rarefaction with the Doubling Rule

**Goal:** Compare diversity across sites at a common sampling completeness, with extrapolation bounded by statistical reliability.

**Approach:** Use iNEXT's coverage-based rarefaction-extrapolation interpolating to the minimum coverage across sites; bound extrapolation at 2x the reference sample size (Chao et al. 2014 *Ecol Monogr* 84:45-67 derive variance bounds that grow superlinearly beyond 2x).

```r
# Type 3 (diversity vs coverage) is the cross-site comparison plot
# nboot=200 minimum for publication-quality CIs (default 50 is too few)
ggiNEXT(result, type = 3) + theme_bw() +
    labs(title = 'Coverage-Based Hill-Number Diversity')

# The doubling rule: endpoint should be at most 2 * max(observed sample size)
# iNEXT default endpoint = 2 * max(sample size); do NOT override above this
# Beyond 2x, the variance estimate grows superlinearly and CIs become unreliable

# To extract numerical results at standardized coverage:
est_95 <- estimateD(abundance_data, datatype = 'abundance',
                    base = 'coverage', level = 0.95)
```

## The Double-Zero Problem and Hellinger Transformation

**Goal:** Compute community dissimilarity in a way that does not treat shared absences as evidence of similarity.

**Approach:** Apply the Hellinger transformation (Legendre & Gallagher 2001 *Oecologia* 129:271-280) before computing Euclidean distance, PCA, or RDA. This is the single most important preprocessing decision for community ordination.

Standard Euclidean distance and Pearson correlation treat "both samples have 0 abundance of species X" as evidence of similarity. In community ecology this is wrong — two deserts both lacking a rainforest species are not thereby similar. Hellinger transformation removes this artifact while preserving total-abundance information.

```r
library(vegan)

# Hellinger transformation: y_ij' = sqrt(y_ij / row_sum_i)
# Euclidean distance on Hellinger-transformed data = Hellinger distance
species_hell <- decostand(community_matrix, method = 'hellinger')

# Now PCA/RDA on the transformed matrix is biologically meaningful
pca_result <- rda(species_hell)

# Alternative: chord transformation (similar effect, slightly different scaling)
species_chord <- decostand(community_matrix, method = 'normalize')
```

## Beta Diversity Partition — Run BOTH Baselga AND Podani

For the broader "multiple meanings of beta diversity" framework, see Anderson et al. 2011 *Ecol Lett* 14:19-28.

**Goal:** Decompose total beta diversity into ecologically interpretable components, acknowledging that the partition is not mathematically unique.

**Approach:** Compute the Baselga partition (Sorensen = turnover + nestedness) with `betapart`, and the alternative Podani/Carvalho partition (richness-difference framework), and report both. The two frameworks give different ecological interpretations of the same data; presenting only one hides ambiguity.

```r
library(betapart)

pa_matrix <- ifelse(community_matrix > 0, 1, 0)

# --- Baselga partition (Baselga 2010 Glob Ecol Biogeogr 19:134-143) ---
# beta.SIM (turnover/Simpson) + beta.SNE (nestedness) = beta.SOR (total Sorensen)
pair_sor <- beta.pair(pa_matrix, index.family = 'sorensen')
multi_sor <- beta.multi(pa_matrix, index.family = 'sorensen')

# --- Abundance-based: Bray-Curtis balanced + gradient ---
# beta.bray.bal (balanced variation; analogous to turnover)
# beta.bray.gra (abundance gradient; analogous to nestedness)
pair_abund <- beta.pair.abund(community_matrix, index.family = 'bray')

# --- Podani/Carvalho framework (different decomposition of same data) ---
# Use betapart::beta.pair with the .fam family option, or carvalho package
# These produce richness-difference components instead of nestedness
# Reporting only Baselga without acknowledging Podani exists is incomplete
```

## Phylogenetic Diversity — Faith's PD, MPD, MNTD with Null-Model Choice

**Goal:** Quantify evolutionary distinctness and phylogenetic community structure with an explicit, justified null model.

**Approach:** Compute Faith's PD on a community matrix and ultrametric phylogeny (Faith 1992 *Biol Conserv* 61:1-10), then SES_MPD and SES_MNTD (Webb 2002 *Annu Rev Ecol Syst* 33:475-505) standardizing against a documented null model. The choice of null model is the dominant scientific decision — different nulls answer different ecological questions.

```r
library(picante)

# Faith's PD: sum of branch lengths spanning the focal species
# include.root=TRUE counts root branch; FALSE excludes (matters for within-clade comparisons)
faith_pd <- pd(community_matrix, phylo_tree, include.root = TRUE)

# SES_MPD with EXPLICIT null model (do not accept default silently)
# null.model='taxa.labels': shuffles species across tree, holds sample richness constant
# null.model='richness': randomizes within sample, preserves species occurrence frequencies
# null.model='independentswap': preserves both row and column sums of community matrix
# Each null answers a different question; document the choice in methods
ses_mpd_taxa <- ses.mpd(community_matrix, cophenetic(phylo_tree),
                        null.model = 'taxa.labels', runs = 999, iterations = 1000)
ses_mpd_indep <- ses.mpd(community_matrix, cophenetic(phylo_tree),
                         null.model = 'independentswap', runs = 999, iterations = 1000)

# SES_MPD < 0 (NRI > 0): phylogenetic clustering
# SES_MPD > 0 (NRI < 0): phylogenetic overdispersion
# DO NOT infer "clustering = environmental filtering" from sign alone;
# Mayfield & Levine 2010 Ecol Lett 13:1085-1093 showed competition can cluster
# when traits track phylogeny and similar species coexist via R*-rule dynamics
```

## Functional Diversity with Trait-Axis Dimensionality Optimization

**Goal:** Compute multidimensional functional diversity without inflating values via correlated trait axes.

**Approach:** Run Maire et al. 2015 *Glob Ecol Biogeogr* 24:728-740 trait-space-quality assessment (mAD metric) to choose the number of PCoA axes that minimizes deviation between original trait distances and axis-reconstructed distances. Use that k for FRic, FEve, FDiv. Default k=2-3 typically biases FRic; optimum is usually 4-6 for typical trait datasets.

```r
library(FD)
library(mFD)

# Trait dissimilarity from a trait matrix (Gower for mixed quantitative/categorical)
trait_dist <- gowdis(trait_matrix)

# mFD optimizes the number of trait axes via Maire 2015 mAD criterion
# Smaller mAD = better fidelity between trait distances and PCoA-axis distances
quality_funct_space <- quality.fspaces(trait_dist, maxdim_pcoa = 10,
                                       fdist_scaling = TRUE, fdendro = NULL)
best_k <- which.min(quality_funct_space$quality_fspaces$mad)
cat('Maire-optimal number of trait axes:', best_k, '\n')

# Compute FD with the optimal k
fd_result <- dbFD(trait_matrix, community_matrix, m = best_k, corr = 'cailliez')
fd_result$FRic   # functional richness (convex hull volume)
fd_result$FEve   # functional evenness
fd_result$FDiv   # functional divergence
```

## Per-Method Failure Modes

### Singleton-driven Chao1 inflation in amplicon data

**Trigger:** Computing Chao1 directly on ASV/OTU tables that include singletons of suspected PCR-error origin.

**Mechanism:** Chao1 assumes singletons are biologically real rare species; under that assumption, more singletons relative to doubletons signals more undiscovered diversity. PCR/sequencing error produces many singletons that look like rare species to Chao1, inflating the bound.

**Symptom:** Chao1 dramatically exceeds observed richness (e.g., Chao1 = 500 from S_obs = 100) with extremely wide confidence intervals; Good's coverage < 0.85 despite > 10,000 reads per sample.

**Fix:** Denoise first (DADA2 / UNOISE3 / swarm v2) to remove PCR-error variants, then either skip Chao1 entirely (Callahan 2017 ASV philosophy) or report observed ASV count + Good's coverage instead. Alternative: report Chao2 from incidence (presence across replicates), which is robust to PCR-error singletons.

### Rarefaction across mismatched assemblage sizes

**Trigger:** Sample-size rarefaction to the smallest assemblage when sites differ in true diversity.

**Mechanism:** A 100-read rarefaction of a 50-species community is nearly saturated (coverage approximately 99%); the same 100 reads from a 500-species community covers approximately 40% of true diversity. The two rarefied values measure different completeness levels.

**Symptom:** Rarefied richness ordering disagrees with intuitive site-by-site comparison; small-sample-size sites appear artificially equal in rarefied diversity to high-diversity sites.

**Fix:** Use coverage-based rarefaction via `estimateD(..., base = 'coverage', level = 0.95)`. The 0.95 coverage target is the modern default; 0.99 for high-precision work, 0.85 minimum for accepting a site into the comparison.

### Hellinger forgotten before PCA/RDA on community data

**Trigger:** Running `rda(species_matrix)` or `prcomp(species_matrix)` on raw species counts without prior transformation.

**Mechanism:** The double-zero problem inflates dissimilarity for sample pairs sharing many absent species. Raw-count PCA projects samples primarily by sequencing depth (PC1 = library size proxy), not by community composition.

**Symptom:** PC1 axis correlates strongly with sample read totals; biological gradients appear only on PC3 or later; ordination is dominated by sites with extreme sample sizes.

**Fix:** `decostand(matrix, method = 'hellinger')` before ordination. Always.

### FRic inflation from collinear trait axes

**Trigger:** Computing FRic from 8-10 raw correlated traits without dimensionality optimization.

**Mechanism:** FRic = convex-hull volume in trait PCoA space. Adding correlated traits increases the apparent dimensionality of the trait space without adding ecological information; convex-hull volume inflates accordingly.

**Symptom:** FRic increases when adding new correlated traits to the analysis; FRic comparisons across studies that use different trait counts are inconsistent.

**Fix:** Run `mFD::quality.fspaces()` and use the mAD-optimal k. Report k in methods.

### SES_MPD interpretation by sign alone

**Trigger:** Reporting "SES_MPD < 0 indicates environmental filtering" without testing alternative explanations.

**Mechanism:** Mayfield & Levine 2010 *Ecol Lett* 13(9):1085-1093 showed competition can produce phylogenetic clustering when ecologically similar species (which tend to be closely related) coexist via R*-rule trait differences. The clustering = filtering / overdispersion = competition dichotomy from Webb 2002 is incomplete.

**Symptom:** Phylogenetic clustering interpretation is challenged at review; trait-similarity vs phylogenetic-similarity correlation has not been examined.

**Fix:** Report SES_MPD with one explicit null model AND test trait conservatism (Blomberg's K, Pagel's lambda) — if traits track phylogeny strongly, clustering may reflect either filtering OR competition; cite Mayfield & Levine 2010 in interpretation.

## Quantitative Thresholds

| Threshold | Value | Source / rationale |
|-----------|-------|-------------------|
| Coverage target for cross-site comparison | C = 0.95 | Postdoc-grade convention; 0.99 for high-precision, 0.85 minimum to accept site |
| iNEXT extrapolation limit | endpoint <= 2x reference sample size | Chao et al. 2014 doubling rule; variance grows superlinearly beyond |
| iNEXT bootstrap floor for CIs | nboot = 200 | Default 50 is too few for publication-quality CIs |
| Chao1 reliability | f2 > 0 AND singletons biological | If f2 = 0, switch to Chao1bc or jackknife1 |
| Hellinger before ordination | Always for community data | Legendre & Gallagher 2001; non-negotiable |
| SES_MPD significance | |SES| > 1.96 (two-tailed) | Standard 0.05 alpha; report alongside explicit null model |
| Functional-diversity axis count k | Maire 2015 mAD-optimal | Default k=2-3 too few; typically 4-6 optimal |
| Phylogenetic-tree requirement | Ultrametric | If not, `ape::chronos()` or BEAST-based dating |

## Common errors

| Error | Cause | Solution |
|-------|-------|----------|
| Chao1 reports infinity or NaN | f2 = 0 (no doubletons) | Use Chao1bc form or jackknife1 |
| iNEXT extrapolation curve flat then explodes | Extrapolated beyond doubling-rule limit | Set endpoint <= 2 * max(sample sizes) |
| Hellinger PCA gives flat results | Decostand not applied or applied to wrong axis | `decostand(matrix, method = 'hellinger')` rows = sites |
| ses.mpd returns all p-values approximately 0.5 | Wrong null model for question (e.g., `taxa.labels` on a single-richness dataset) | Choose null that varies the quantity tested |
| FRic dramatically increases with new traits | Correlated trait axes inflating convex hull | Optimize axis count with mFD::quality.fspaces |
| beta.multi returns NaN for nestedness | Communities entirely disjoint (no shared species) | beta.SNE -> 0 trivially; interpret as pure turnover |
| Bray-Curtis flagged for triangle-inequality violation | Bray-Curtis is not a true metric | Switch to Sorensen (a metric) for downstream methods requiring metricity |

## References

- Chao A, Jost L (2012) Coverage-based rarefaction and extrapolation. *Ecology* 93(12):2533-2547. doi:10.1890/11-1952.1
- Chao A, Gotelli NJ, Hsieh TC, Sander EL, Ma KH, Colwell RK, Ellison AM (2014) Rarefaction and extrapolation with Hill numbers. *Ecol Monogr* 84(1):45-67. doi:10.1890/13-0133.1
- Hsieh TC, Ma KH, Chao A (2016) iNEXT: an R package for rarefaction and extrapolation of species diversity. *Methods Ecol Evol* 7(12):1451-1456. doi:10.1111/2041-210X.12613
- Jost L (2006) Entropy and diversity. *Oikos* 113(2):363-375. doi:10.1111/j.2006.0030-1299.14714.x
- Faith DP (1992) Conservation evaluation and phylogenetic diversity. *Biol Conserv* 61(1):1-10. doi:10.1016/0006-3207(92)91201-3
- Legendre P, Gallagher ED (2001) Ecologically meaningful transformations for ordination. *Oecologia* 129(2):271-280. doi:10.1007/s004420100716
- Baselga A (2010) Partitioning the turnover and nestedness components of beta diversity. *Glob Ecol Biogeogr* 19(1):134-143. doi:10.1111/j.1466-8238.2009.00490.x
- Anderson MJ, Crist TO, Chase JM et al. (2011) Navigating the multiple meanings of beta diversity. *Ecol Lett* 14(1):19-28. doi:10.1111/j.1461-0248.2010.01552.x
- Webb CO, Ackerly DD, McPeek MA, Donoghue MJ (2002) Phylogenies and community ecology. *Annu Rev Ecol Syst* 33:475-505. doi:10.1146/annurev.ecolsys.33.010802.150448
- Mayfield MM, Levine JM (2010) Opposing effects of competitive exclusion on phylogenetic community structure. *Ecol Lett* 13(9):1085-1093. doi:10.1111/j.1461-0248.2010.01509.x
- Maire E, Grenouillet G, Brosse S, Villeger S (2015) How many dimensions are needed to accurately assess functional diversity? *Glob Ecol Biogeogr* 24(6):728-740. doi:10.1111/geb.12299
- Chao A (1984) Nonparametric estimation of the number of classes in a population. *Scand J Stat* 11(4):265-270
- Gloor GB, Macklaim JM, Pawlowsky-Glahn V, Egozcue JJ (2017) Microbiome datasets are compositional: and this is not optional. *Front Microbiol* 8:2224. doi:10.3389/fmicb.2017.02224

## Related Skills

- ecological-genomics/edna-metabarcoding - Generate ASV/species tables prior to diversity analysis
- ecological-genomics/community-ecology - Constrained ordination, indicator species, PERMANOVA on transformed data
- microbiome/diversity-analysis - 16S clinical microbiome diversity metrics with compositional considerations
- data-visualization/ggplot2-fundamentals - Customize diversity plots and rarefaction curves
- phylogenetics/tree-io - Ultrametric tree preparation for PD/MPD/MNTD
