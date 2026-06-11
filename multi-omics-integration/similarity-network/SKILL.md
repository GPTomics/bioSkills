---
name: bio-multi-omics-similarity-network
description: Stratifies patients into multi-omics subtypes by building one patient-by-patient similarity network per omic, fusing them with SNF's cross-network diffusion, and spectral-clustering the fused graph - then defending the clusters with stability, survival separation, and replication. Covers why spectral clustering always returns the requested cluster count so a subtype is a claim not a discovery, why the eigengap is a graph property not a biological truth, why fusion is not automatically better than the best single omic, why SNF needs complete data while NEMO handles mosaic cohorts, and the SNFtool API gotchas (dist2 returns squared distance, affinityMatrix width is sigma, spectralClustering K is the cluster count). Use when discovering patient subtypes from multiple omics, choosing a cluster number, validating subtypes, or handling partial multi-omic data. For feature-space factors see mofa-integration; for supervised signatures see mixomics-analysis; for survival see clinical-biostatistics/survival-analysis.
tool_type: r
primary_tool: SNFtool
---

## Version Compatibility

Reference examples tested with: SNFtool 2.3+, igraph 2.0+, pheatmap 1.0+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('SNFtool')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

SNFtool argument names are easy to misread: `affinityMatrix`'s width argument is `sigma` (not alpha or mu), `dist2` returns SQUARED Euclidean distance (take the square root), and `spectralClustering`'s `K` is the number of CLUSTERS, a name collision with the K-neighbors argument of `affinityMatrix`/`SNF`.

# Similarity Network Fusion

**"Stratify my patients using multi-omics data"** -> Fuse per-omic patient-similarity networks into one graph and spectral-cluster it - because the clustering always returns the number of subtypes requested, so the subtype is a claim to defend, not a discovery.
- R: `SNF()` to fuse networks, `spectralClustering()` to partition, then validate

Scope: patient-similarity-space integration - per-omic affinity networks, SNF fusion, spectral clustering into candidate subtypes, cluster-number choice, and the stability/survival/replication defense, plus the integrative-clustering landscape (NEMO, PINS, iCluster, CIMLR, intNMF, consensus). Feature-space latent factors -> mofa-integration. Supervised feature signatures -> mixomics-analysis. Survival mechanics -> clinical-biostatistics/survival-analysis. Single-cell graph clustering -> single-cell/clustering.

## The Single Most Important Modern Insight -- Spectral Clustering Always Returns the Requested Number of Clusters, So a Subtype Is a Claim to Defend, Not a Thing Discovered

Ask the spectral step for four clusters and it returns four, whether or not four real patient groups exist; the eigengap suggesting "four" means only that the fused graph has a roughly four-block shape, not that the disease has four subtypes. SNF's fusion is genuinely powerful - its cross-network diffusion reinforces patient pairs that multiple omics agree on and erodes omic-specific noise - but power to draw a boundary is not evidence the boundary is real. Three defenses, none of which the algorithm supplies:

1. **Stability.** Resample the patients, re-cluster, and measure agreement (consensus index / NMI across subsamples). A real subtype recurs; an artifact does not. Report stability across a fixed (K, sigma) grid rather than tuning those hyperparameters to a target.
2. **Survival separation after adjustment.** Subtypes should differ in outcome in a Cox model adjusted for known prognostic covariates (stage, age, grade), so a "subtype" is not just re-encoded stage. The discovery-cohort p is supporting evidence, not proof - and a permutation p is safer than the chi-square approximation at small arm sizes.
3. **Replication and the fused-versus-best-single-omic check.** Multi-omics does not consistently beat the best single omic (Rappoport and Shamir 2018), and SNF's diffusion can dilute a signal carried by one omic. Benchmark the fused clustering against each single-omic clustering, and replicate the subtypes in an independent cohort.

The honest report names the hyperparameters, shows the clusters are not knife-edge-sensitive to them, and treats the discovery survival p as supporting evidence, not the finding.

## Tool Taxonomy (Integrative Clustering)

| Tool | Citation | Mechanism | When |
|------|----------|-----------|------|
| SNF (SNFtool) | Wang 2014 *Nat Methods* 11:333 | per-omic affinity + cross-network diffusion -> spectral | complete data; non-linear consensus reinforcing multi-omic agreement; the default |
| NEMO | Rappoport 2019 *Bioinformatics* 35:3348 | per-omic relative similarity, AVERAGED (no iteration) -> spectral | PARTIAL/mosaic data; fast; comparable accuracy on full data |
| PINS / PINSPlus | Nguyen 2017 *Genome Res* 27:2025; Nguyen 2019 *Bioinformatics* 35:2843 | perturbation clustering; keeps partitions robust to noise; auto-picks C | when partition stability is the priority |
| iCluster / iClusterPlus | Shen 2009 *Bioinformatics* 25:2906; Mo 2013 *PNAS* 110:4245 | model-based joint latent-variable + feature selection | want a generative model and feature selection; tends to pick few clusters |
| CIMLR | Ramazzotti 2018 *Nat Commun* 9:4453 | multiple-kernel learning per omic -> k-means | one Gaussian kernel per omic too rigid; strong survival results |
| intNMF | Chalise 2017 *PLoS One* 12:e0176278 | joint non-negative matrix factorization | non-negative data; parts-based factorization + clustering |
| Consensus / COCA | Monti 2003 *Mach Learn* 52:91; Hoadley 2014 *Cell* 158:929 | cluster each omic, then cluster the matrix-of-clusters | late integration; want each omic's clustering visible and a consensus |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Complete multi-omics, want a non-linear consensus partition | SNF + spectralClustering | cross-network diffusion reinforces multi-omic agreement |
| Some patients missing an omic (mosaic/partial) | NEMO | built for partial data; no imputation, no patient loss |
| Partition robustness/auto cluster number is the priority | PINSPlus | perturbation clustering builds stability in |
| Want a generative model and integrated feature selection | iCluster / iClusterBayes | model-based joint latent variable |
| Need feature-level interpretation (which genes define a subtype) | -> mofa-integration / mixomics-analysis | SNF has no feature model; feature-space tools live there |
| Validate subtypes against outcome | -> clinical-biostatistics/survival-analysis | Cox / log-rank / KM mechanics |
| SNF underperforms on survival in a benchmark | MCCA (survival) or rMKL-LPP (clinical enrichment) | topped those criteria in Rappoport and Shamir 2018 |
| Which method at all / paired vs mosaic | -> integration-design | the correspondence and method decision |

## Build the Per-Omic Affinity Networks

**Goal:** Turn each omic into a patient-by-patient similarity network on a common scale, collapsing the high-dimensional feature space into an n-by-n object so feature count buys no votes.

**Approach:** Standardize each continuous omic per feature, compute the (square-rooted) Euclidean distance, then apply the local-scaled Gaussian kernel. SNF requires every patient to have every omic, so intersect to common samples first and report how many that drops.

```r
library(SNFtool)

K <- 20        # neighbors defining the local kernel bandwidth; SNFtool guidance 10-30; changes cluster count
sigma <- 0.5   # kernel width multiplier (affinityMatrix's third arg, named sigma not alpha); guidance 0.3-0.8
t_iter <- 20   # cross-diffusion iterations; converges by ~10-20

norm_views <- lapply(list(rna=rna, meth=meth, mirna=mirna), standardNormalization)   # per-feature z-score before distance
dists <- lapply(norm_views, function(x) dist2(x, x)^(1/2))                            # dist2 returns SQUARED distance
affinities <- lapply(dists, function(d) affinityMatrix(d, K, sigma))
```

## Fuse and Cluster

**Goal:** Fuse the per-omic networks into one graph and partition it, treating the cluster number as the central claim rather than a nuisance parameter.

**Approach:** Run SNF's cross-diffusion, read the four cluster-number estimates the package returns (not as truth but as plausibility), then spectral-cluster. The package itself warns the estimates cannot guarantee accuracy.

```r
fused <- SNF(affinities, K, t_iter)
estimateNumberOfClustersGivenGraph(fused, NUMC=2:8)   # returns FOUR estimates: K1/K12 (eigengap), K2/K22 (rotation cost)
clusters <- spectralClustering(fused, K=4, type=3)    # here K is the CLUSTER COUNT (not neighbors); type 3 = Ng-Jordan-Weiss default
```

## Defend the Subtypes

**Goal:** Show the clusters are stable, separate outcome, and are not just the best single omic before calling them subtypes.

**Approach:** Compare the fused clustering against each single-omic clustering, assess stability under resampling, and rank the post-hoc feature attribution with the package function rather than a hand-rolled test. Survival mechanics are routed out.

```r
concordanceNetworkNMI(c(affinities, list(fused)), C=4)   # NMI among per-omic and fused clusterings: did fusion beat the best single omic?
feat_rank <- rankFeaturesByNMI(norm_views, fused)        # POST-HOC attribution: features that track the clusters, not a model that made them
```

Validate survival separation in a covariate-adjusted Cox model and report events per arm (clinical-biostatistics/survival-analysis owns the mechanics); assess stability by resampling patients and re-clustering, reporting agreement across subsamples. To assign a new patient to an existing subtype without re-clustering, use `groupPredict(train_views, test_views, groups, K=20, method=1)` (label propagation). For a mosaic cohort, switch to NEMO rather than dropping the incomplete patients.

## Per-Method Failure Modes

### Subtype count reported as a discovery
**Trigger:** "I found 4 subtypes" with no defense. **Mechanism:** spectral clustering returns the requested C regardless of structure. **Symptom:** a clean-looking partition that does not replicate. **Fix:** require eigengap plausibility, resampling stability, covariate-adjusted survival, and external replication before claiming a subtype count.

### Hyperparameters tuned to a label, then claimed
**Trigger:** grid-searching K and sigma to maximize NMI against known labels. **Mechanism:** in real discovery there are no labels, so tuning to NMI is circular. **Symptom:** a result that only holds at the chosen (K, sigma). **Fix:** fix or pre-register a small (K, sigma) grid and show the clustering is stable across it; report sensitivity as a result.

### Complete-data requirement ignored
**Trigger:** `Reduce(intersect, ...)` silently dropping patients missing an omic. **Mechanism:** SNF's cross-diffusion multiplies aligned n-by-n matrices, so it needs complete data. **Symptom:** a decimated, biased cohort. **Fix:** report the dropped count and bias; use NEMO for partial data; do not impute a whole omic to keep a patient.

### Fusion assumed to beat the best single omic
**Trigger:** reporting the fused clustering without a single-omic comparison. **Mechanism:** fusion's diffusion can dilute a signal carried by one omic; multi-omics is not consistently better (Rappoport and Shamir). **Symptom:** a fused result no better than the best layer. **Fix:** benchmark fused vs each single omic with `concordanceNetworkNMI` and per-omic survival.

### Survival p taken at face value
**Trigger:** a log-rank p on discovery clusters. **Mechanism:** C, K, sigma were chosen partly to get separable groups, and arms are small. **Symptom:** an over-optimistic p that does not replicate. **Fix:** adjust for prognostic covariates, report events per arm, use a permutation p, and require replication.

### Feature attribution treated as the model
**Trigger:** presenting ranked features as what generated the clusters. **Mechanism:** SNF has no feature-level model; attribution is post-hoc. **Symptom:** causal claims SNF cannot support. **Fix:** use `rankFeaturesByNMI` and present it as post-hoc characterization; for a feature model use MOFA/DIABLO.

### SNFtool API misread
**Trigger:** treating `dist2` output as Euclidean, passing alpha to `affinityMatrix`, or reading `spectralClustering`'s K as neighbors. **Mechanism:** `dist2` is squared, the width arg is `sigma`, and that `K` is the cluster count. **Symptom:** distorted affinities or the wrong number of clusters. **Fix:** take `dist2(...)^(1/2)`, pass `sigma`, and read the K name in context.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| K neighbors ~20 (range 10-30) | Wang 2014 *Nat Methods* 11:333; SNFtool docs | sets the local kernel bandwidth; small K fragments, large K over-smooths and changes cluster count |
| sigma ~0.5 (range 0.3-0.8) | SNFtool docs | kernel width multiplier; wider blurs clusters, narrower sharpens noise |
| t iterations ~10-20 | Wang 2014 *Nat Methods* 11:333 | cross-diffusion converges; more iterations do little past convergence |
| `estimateNumberOfClustersGivenGraph` returns FOUR estimates | SNFtool docs | eigengap (K1/K12) and rotation cost (K2/K22); plausibility not proof |
| Fused must beat the best single omic to justify fusion | Rappoport and Shamir 2018 *Nucleic Acids Res* 46:10546 | multi-omics is not consistently better; check explicitly |
| Stability across a fixed (K, sigma) grid + permutation survival p | Monti 2003 *Mach Learn* 52:91 | resampling stability and a permutation p guard against artifact subtypes |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Distorted affinities / wrong-scale distances | `dist2` output used as Euclidean | take `dist2(...)^(1/2)` |
| Clusters change unexpectedly | `affinityMatrix` width passed as alpha or wrong arg | the third argument is `sigma` |
| Wrong number of clusters | `spectralClustering`'s K read as neighbors | that `K` is the cluster count |
| Function `dist2` not found after first call | a variable named `dist2` shadowing the function | name distance variables differently (d1, d2) |
| Many patients silently dropped | complete-data intersection on a mosaic cohort | report the drop; use NEMO for partial data |
| Hand-rolled feature ranking | reimplementing attribution with `aov` | use `rankFeaturesByNMI(list_of_views, fused)` |

## References

- Wang B, Mezlini AM, Demir F, et al. 2014. Similarity network fusion for aggregating data types on a genomic scale. *Nat Methods* 11:333-337.
- Rappoport N, Shamir R. 2018. Multi-omic and multi-view clustering algorithms: review and cancer benchmark. *Nucleic Acids Res* 46:10546-10562.
- Rappoport N, Shamir R. 2019. NEMO: cancer subtyping by integration of partial multi-omic data. *Bioinformatics* 35:3348-3356.
- Nguyen T, Tagett R, Diaz D, Draghici S. 2017. A novel approach for data integration and disease subtyping. *Genome Res* 27:2025-2039.
- Nguyen H, Shrestha S, Draghici S, Nguyen T. 2019. PINSPlus: a tool for tumor subtype discovery in integrated genomic data. *Bioinformatics* 35:2843-2846.
- Shen R, Olshen AB, Ladanyi M. 2009. Integrative clustering of multiple genomic data types using a joint latent variable model with application to breast and lung cancer subtype analysis. *Bioinformatics* 25:2906-2912.
- Mo Q, Wang S, Seshan VE, et al. 2013. Pattern discovery and cancer gene identification in integrated cancer genomic data. *PNAS* 110:4245-4250.
- Ramazzotti D, Lal A, Wang B, Batzoglou S, Sidow A. 2018. Multi-omic tumor data reveal diversity of molecular mechanisms that correlate with survival. *Nat Commun* 9:4453.
- Chalise P, Fridley BL. 2017. Integrative clustering of multi-level 'omic data based on non-negative matrix factorization algorithm. *PLoS One* 12:e0176278.
- Monti S, Tamayo P, Mesirov J, Golub T. 2003. Consensus clustering: a resampling-based method for class discovery and visualization of gene expression microarray data. *Mach Learn* 52:91-118.
- von Luxburg U. 2007. A tutorial on spectral clustering. *Stat Comput* 17:395-416.

## Related Skills

- integration-design - The method-selection and paired-vs-mosaic decision
- mofa-integration - Feature-space latent factors (which molecular axes matter)
- mixomics-analysis - Supervised feature signatures for known classes
- data-harmonization - Per-omic scaling before building distances
- clinical-biostatistics/survival-analysis - Survival validation of discovered subtypes
- single-cell/clustering - Graph clustering of cells (different object)
- pathway-analysis/go-enrichment - Enrichment of subtype-associated features
- workflows/multi-omics-pipeline - End-to-end multi-omics integration pipeline
