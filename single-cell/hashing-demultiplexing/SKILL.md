---
name: bio-single-cell-hashing-demultiplexing
description: Assign cells to their sample of origin from cell or nucleus hashing (CITE-seq HTOs, MULTI-seq lipid/cholesterol tags, CellPlex CMOs) and call cross-sample doublets using Seurat HTODemux/MULTIseqDemux, hashsolo, demuxEM, GMM-Demux, and demuxmix. Use when assigning pooled hashed cells back to their sample, calling cross-sample doublets from HTO counts, choosing a demultiplexing method, deciding between hashtag and genetic demultiplexing, or rescuing an oversized Negative pile from weak HTO staining or ambient spillover.
tool_type: mixed
primary_tool: Seurat
---

## Version Compatibility

Reference examples tested with: Seurat 5.0+, scanpy 1.10+, pegasus 1.8+, demuxmix 1.4+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Hashtag Demultiplexing and Cross-Sample Doublet Calling

**"Which sample did each cell come from, and which barcodes are cross-sample doublets?"** -> Classify every cell's HTO count vector to the one tag that dominates its background, and flag cells where two tags are both high.
- R: `Seurat::HTODemux` (antibody HTOs), `Seurat::MULTIseqDemux` (MULTI-seq lipid tags), `demuxmix` (regression mixture, robust to bad staining)
- Python: `scanpy.external.pp.hashsolo` (Bayesian), `pegasus.demultiplex` / demuxEM (background from empty droplets)
- CLI: `GMM-Demux` (Gaussian mixture with explicit multiplet accounting)

## Governing principle

In cell or nucleus hashing, each sample is labeled before pooling with a unique oligo-tagged reagent: a barcoded antibody (CITE-seq HTO), a lipid- or cholesterol-modified oligo (MULTI-seq), or a CellPlex CMO. A true singlet's HTO count vector is dominated by ONE tag standing well above a background of ambient and spillover counts, so sample assignment reduces to one question per cell: which tag, if any, exceeds that cell's background. Cross-sample doublets fall out directly - two tags both high - which is the decisive advantage of hashing over expression-only doublet detection.

Hashtag demux, genetic demux, and expression-doublet detection answer three DIFFERENT questions and should be combined, not substituted. Hashtag and genetic demux both assign samples and catch cross-sample doublets, but expression-doublet detection (single-cell/doublet-detection) catches within-sample and homotypic doublets that hashing and genetics are blind to, because two cells from the same sample carry the same tag and the same genotype. Conversely, expression methods miss cross-sample doublets when the two samples are transcriptionally similar. The cross-sample doublet rate also calibrates the expected TOTAL doublet rate: with two pooled samples within- and cross-sample doublets are equally frequent, while with k samples cross-sample doublets dominate and within-sample ones fall to about 1/k of all doublets, so hashing still misses that within-sample fraction (expression-doublet detection stays necessary) and a hashing doublet rate far below the expression-doublet rate is a red flag that staining or thresholds are off.

The hard part is the background, not the dominant tag. Ambient HTO from lysed cells, spillover between tags, staining failure, and batch differences in tag-capture efficiency all inflate the background and grow the "Negative" pile (real cells whose true tag never cleared background, distinct from true empty droplets removed upstream). Nucleus hashing is harder than whole-cell because tag capture is lower. Methods differ mainly in how they model that background.

## Choosing a demultiplexing modality

| Modality | Needs | Cross-sample doublets | Cannot do | Tools |
|----------|-------|-----------------------|-----------|-------|
| Hashtag/HTO (this skill) | HTO/lipid/CMO library at pooling | Yes (two tags high) | Nothing without a hashing library; sensitive to staining and ambient | HTODemux, MULTIseqDemux, hashsolo, demuxEM, GMM-Demux, demuxmix |
| Genetic (natural SNPs) | >=2 distinct genotypes, no hashing | Yes | Cannot separate same-donor samples (identical genotype) | demuxlet, freemuxlet, souporcell, vireo |
| Expression doublet (orthogonal) | Just the GEX matrix | No (misses when samples similar) | Catches within-sample/homotypic doublets the other two miss | scDblFinder, Scrublet (single-cell/doublet-detection) |

Genetic demux is the fallback when no hashing was done but samples come from different donors; it cannot resolve multiple samples from one donor, which hashing can. Pair whichever sample-assignment method applies with expression-doublet detection for the doublets it cannot see.

## Choosing a hashtag caller

| Method | Model | Use when | Fails when |
|--------|-------|----------|------------|
| HTODemux (Seurat) | k-medoids cluster per HTO + negative-distribution quantile | Standard antibody HTO, clean bimodal staining, Seurat workflow | Weak/low-depth staining or heavy ambient; clustering unstable on near-zero HTOs |
| MULTIseqDemux (Seurat) | Per-HTO KDE, threshold between maxima, quantile sweep | MULTI-seq lipid/cholesterol tags; want autoThresh to optimize the quantile | Few cells; unimodal density when one tag dominates |
| hashsolo (scanpy/solo) | Bayesian over negative/singlet/doublet | Few hashtags (works at 2), many negatives, scanpy-native pipeline | Very low signal; priors mis-set for the actual doublet rate |
| demuxEM (pegasus) | EM with background estimated from empty droplets | High ambient; nucleus hashing; raw matrix with empties available | Empty droplets filtered out before calling; very sparse signal |
| GMM-Demux (CLI) | Gaussian mixture on normalized HTO, explicit MSM multiplets | Want explicit multiplet accounting or experiment planning | Non-Gaussian background; poor per-tag separation |
| demuxmix (R) | Negative-binomial regression mixture, optional RNA covariate | Bad/variable staining, batch tag-efficiency differences | Very few cells per mixture component |

The EM and regression methods (demuxEM, demuxmix) model the background explicitly and are the robust choice when staining is marginal, and they handle a two-tag pool as readily as a many-tag one (their failure mode is too few cells per mixture component, not too few tags), so weak staining even at two tags routes to demuxmix or demuxEM rather than hashsolo; HTODemux and MULTIseqDemux are fast defaults for clean data; hashsolo handles few hashes and many negatives. When callers disagree, run a consensus (cellhashR wraps several callers) and verify current best practice against installed docs before trusting any single call.

## Normalize HTO counts before calling

**Goal:** Put HTO counts on a scale where the dominant tag separates from background.

**Approach:** Apply the centered log-ratio (CLR) transform to the HTO assay. CLR margin=1 normalizes the tags within each cell and is the Seurat default that the canonical HTO vignette uses; margin=2 normalizes each tag across cells and is a common alternative for HTO/ADT because it corrects per-tag capture-efficiency differences. Choose deliberately and compare both rather than blindly accepting the default.

```r
library(Seurat)

hto <- CreateSeuratObject(counts = gex_counts)
hto[['HTO']] <- CreateAssay5Object(counts = hto_counts)
hto <- NormalizeData(hto, assay = 'HTO', normalization.method = 'CLR', margin = 2)
```

## Classify samples and doublets with HTODemux (R)

**Goal:** Assign each cell to a single HTO or label it a cross-sample doublet or Negative.

**Approach:** Cluster cells per HTO, model the low-count (negative) cluster, and call a cell positive for any tag whose count exceeds the `positive.quantile` of that negative distribution; one positive is a singlet, two or more a doublet, none a Negative.

```r
hto <- HTODemux(hto, assay = 'HTO', positive.quantile = 0.99)

table(hto$HTO_classification.global)        # Singlet / Doublet / Negative
table(hto$hash.ID)                          # per-sample singlet counts + Doublet + Negative
singlets <- subset(hto, subset = HTO_classification.global == 'Singlet')
```

`positive.quantile = 0.99` is the quantile of the inferred negative distribution above which a cell counts as positive; raise it to be stricter (fewer false singlets, more Negatives), lower it to rescue cells when staining is weak. `HTO_classification.global` holds Singlet/Doublet/Negative; `hash.ID` holds the sample name (or Doublet/Negative) and becomes the active identity.

## Demultiplex MULTI-seq tags with MULTIseqDemux (R)

**Goal:** Classify MULTI-seq lipid/cholesterol-tagged samples, optimizing the threshold automatically.

**Approach:** For each tag, find the threshold between the two density maxima; with `autoThresh=TRUE`, sweep the quantile over `qrange` to maximize the number of singlets.

```r
hto <- MULTIseqDemux(hto, assay = 'HTO', autoThresh = TRUE)
table(hto$MULTI_ID)                          # sample / Doublet / Negative
```

`MULTI_ID` carries the per-cell call. Use a fixed `quantile = 0.7` instead of `autoThresh` only when the automated sweep over- or under-calls on a particular dataset.

## Demultiplex in scanpy with hashsolo (Python)

**Goal:** Bayesian sample assignment that behaves with few hashtags and many negatives.

**Approach:** Place raw HTO counts as columns in `adata.obs`, then run hashsolo with priors over the negative, singlet, and doublet hypotheses; the doublet prior should track the expected loading doublet rate.

```python
import scanpy as sc
import scanpy.external as sce

hto_cols = ['HTO_A', 'HTO_B', 'HTO_C', 'HTO_D']
adata.obs[hto_cols] = hto_counts_df[hto_cols]
sce.pp.hashsolo(adata, cell_hashing_columns=hto_cols, priors=(0.01, 0.8, 0.19))

adata.obs['Classification'].value_counts()   # barcode name / 'Negative' / 'Doublet'
singlets = adata[~adata.obs['Classification'].isin(['Negative', 'Doublet'])].copy()
```

`priors` are ordered [negative, singlet, doublet]; raise the doublet prior for higher loading. Output columns include `Classification`, `most_likely_hypothesis`, and the per-hypothesis probabilities.

## Model ambient background explicitly (Python / R)

**Goal:** Recover correct calls when ambient HTO or weak staining inflates the background.

**Approach:** demuxEM estimates the background from empty droplets before assigning signal; demuxmix fits a negative-binomial regression mixture using the number of detected genes as a covariate, both of which are more robust than a fixed quantile.

```python
import pegasus as pg

pg.estimate_background_probs(hashing_data)
pg.demultiplex(rna_data, hashing_data, min_signal=10.0)
rna_data.obs['demux_type'].value_counts()    # singlet / doublet / unknown
rna_data.obs['assignment']                    # sample name per cell
```

```r
library(demuxmix)

dmm <- demuxmix(as.matrix(hto_counts), rna = num_detected_genes)
calls <- dmmClassify(dmm)                      # HTO assignment + Type (singlet/multiplet/negative/uncertain)
```

`min_signal=10.0` marks cells with too little signal as unknown; lower it to rescue low-capture nucleus hashing, raise it for cleaner singlets. demuxmix's RNA covariate is what makes it robust to per-tag staining differences.

## Threshold and parameter reference

| Parameter | Default | Rationale and when to change |
|-----------|---------|------------------------------|
| HTODemux positive.quantile | 0.99 | Quantile of the negative distribution defining "positive"; raise for stricter calls (more Negatives), lower to rescue weak staining |
| NormalizeData CLR margin | 1 (Seurat default); 2 common for HTO | margin=2 normalizes each tag across cells, correcting per-tag capture bias; pick per the staining and verify |
| MULTIseqDemux quantile / autoThresh | 0.7 / FALSE | autoThresh sweeps the quantile to maximize singlets; use when a fixed threshold over- or under-calls |
| hashsolo priors | (0.01, 0.8, 0.19) | [negative, singlet, doublet]; the doublet prior should track expected loading doublets (~0.8% per 1000 cells on 10x) |
| demuxEM min_signal | 10.0 | Cells below this signal are unknown; lower for low-capture nuclei, raise for cleaner singlets |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Huge Negative pile, few singlets | Weak staining or high ambient inflating background; quantile too strict | Lower positive.quantile / min_signal; switch to demuxEM (empty-droplet background) or demuxmix (RNA covariate) |
| Cross-sample doublet rate near zero but expression doublets high | Hashing thresholds too loose, or doublet prior too low | Tighten the quantile; raise hashsolo doublet prior; reconcile against the expected loading doublet rate |
| HTODemux errors on a zero-count cluster | Cells with all-zero HTO counts cluster together | Filter cells with no HTO counts before HTODemux; check the HTO matrix barcodes match the GEX cells |
| Calls flip with normalization choice | CLR margin=1 vs margin=2 shifts per-tag thresholds | Choose margin deliberately (margin=2 corrects tag-efficiency bias); compare both and inspect ridge plots |
| Genetic demux cannot split two samples | Both samples are the same donor (identical genotype) | Use hashtag demux; genetic methods cannot separate same-donor samples |
| Cross-sample doublets present but homotypic doublets remain | Hashing is blind to within-sample doublets | Run expression-doublet detection (single-cell/doublet-detection) in addition |
| Nucleus hashing yields mostly Negatives | Lower tag capture in nuclei than whole cells | Use demuxEM (designed for nuclei); lower min_signal; expect a larger Negative fraction |
| Looks clean (low Negatives, low doublets) but nearly all cells are one sample | Staining failure where one tag dominates all cells (mispipetted/over-concentrated antibody, or all samples got the same tag) | Sanity-check the per-tag singlet distribution against the expected pooling; one tag capturing nearly all cells means staining failed even though Negatives look low |
| One sample silently lost or contaminating while others demultiplex fine | A single antibody failed to stain, so its cells fall into Negative or misassign to the nearest-ambient tag | Check each tag has a non-trivial positive population; one near-zero tag means a failed antibody dropped or misassigned that sample |
| The rare sample in unequal pooling is under-recovered | Very unequal pooling (e.g. 80/10/10) leaves the minority tag too few positives to form a clean cluster or negative distribution | Inspect per-tag ridge plots; consider demuxmix for the minority tag, whose regression mixture is more stable on small components |
| Many cells flagged generic "Doublet" | Cells positive for 3+ tags collapsed to one label, hiding over-loading or heavy ambient | Inspect the multiplet tag-count distribution (GMM-Demux MSM); 3+ tags high is a run-quality diagnostic, not an ordinary 2-cell doublet |

## Related Skills

- single-cell/doublet-detection - Expression-based within-sample doublet calling that complements cross-sample hashing doublets
- single-cell/preprocessing - Filter empty droplets and QC the cells before and after demultiplexing
- single-cell/batch-integration - Integrate the demultiplexed per-sample data; covers genetic demultiplexing as an alternative
- single-cell/multimodal-integration - HTOs are an ADT-like modality; the CLR normalization here parallels CITE-seq ADT handling
- single-cell/clustering - Cluster the recovered singlets after sample assignment

## References

- Stoeckius et al. 2018, Genome Biol 19:224 - Cell Hashing; barcoded antibodies for multiplexing and doublet detection.
- McGinnis et al. 2019, Nat Methods 16:619-626 - MULTI-seq; lipid- and cholesterol-tagged-oligo sample multiplexing.
- Kang et al. 2018, Nat Biotechnol 36:89-94 - demuxlet; genetic demultiplexing from natural variation.
- Heaton et al. 2020, Nat Methods 17:615-620 - souporcell; genotype clustering without reference genotypes.
- Huang et al. 2019, Genome Biol 20:273 - vireo; Bayesian genetic demultiplexing without a genotype reference.
- Bernstein et al. 2020, Cell Syst 11(1):95-101 - Solo and hashsolo; Bayesian hashing demultiplexing.
- Gaublomme et al. 2019, Nat Commun 10:2907 - demuxEM; nuclei multiplexing with background estimated from empty droplets.
- Xin et al. 2020, Genome Biol 21:188 - GMM-Demux; Gaussian mixture with multi-sample-multiplet accounting.
- Klein 2023, Bioinformatics 39(8):btad481 - demuxmix; negative-binomial regression mixture robust to staining differences.
