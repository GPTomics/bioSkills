---
name: bio-spatial-transcriptomics-spatial-deconvolution
description: Estimates per-spot cell type composition of spatial transcriptomics mixtures (Visium, Slide-seq, Stereo-seq) from an scRNA-seq reference with cell2location, RCTD, SPOTlight, stereoscope, SpatialDWLS, or reference-free STdeconvolve. Use when deciding whether a platform even needs deconvolution (the resolution fork -- a 55um Visium spot is a 1-10-cell MIXTURE -> deconvolve, but a Xenium/MERFISH/CosMx cell is already single -> segment instead, and running deconvolution there invents fractions that do not exist); choosing cell2location (absolute abundance) vs RCTD/SPOTlight/stereoscope/SpatialDWLS (proportions only) by output and runtime; matching the scRNA reference to tissue and condition (the reference IS the result -- a missing cell type is silently misassigned to its nearest neighbor with no error flag); and handling compositional outputs that sum to 1 with CLR/ILR rather than naive per-type t-tests.
tool_type: python
primary_tool: cell2location
---

## Version Compatibility

Reference examples tested with: cell2location 0.1.4+, scvi-tools 1.0+, scanpy 1.10+, anndata 0.10+, numpy 1.26+

RCTD/SPOTlight/STdeconvolve are R packages (spacexr, SPOTlight, STdeconvolve via Bioconductor); call them from R or via rpy2.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Spatial Deconvolution

**"What cell types are in each spot?"** -> Decompose a multi-cell capture spot into the fractions (or absolute numbers) of each cell type, using an annotated scRNA-seq reference as the basis.
- Python: cell2location (`RegressionModel` -> `Cell2location`), stereoscope/DestVI (scvi-tools), Tangram (`tg.map_cells_to_space`), STdeconvolve (reference-free, R)
- R: RCTD (`spacexr::create.RCTD`/`run.RCTD`), SPOTlight, SpatialDWLS (Giotto), CARD

The operation BRANCHES on the platform before any tool is chosen. The first question is not "which method" but "does this data even contain mixtures?" -- answered by the resolution fork below.

## The resolution fork

Deconvolution recovers the cell-type composition of a MIXTURE. Whether a mixture exists is set by the capture unit size relative to a mammalian cell (~8-30um diameter), and it sorts every dataset into three regimes. Choosing the wrong regime is the most expensive error in spatial analysis -- far more damaging than picking the second-best method within a regime.

| Platform | Unit size | Single cell? | Regime |
|----------|-----------|--------------|--------|
| Visium v1/v2 | 55um spot, 100um pitch | No (~1-10+ cells/spot) | DECONVOLVE |
| GeoMx DSP | region of interest | No (many cells) | DECONVOLVE (SpatialDecon) |
| Visium HD | 2um bins, analyzed at 8/16um | 8um still mixes cells | AMBIGUOUS (reconstruct OR deconvolve) |
| Stereo-seq | ~220nm spots, binned (bin50/bin100) | binned to cell scale | AMBIGUOUS |
| Slide-seqV2 | 10um beads | near single-cell | AMBIGUOUS (RCTD doublet-mode common) |
| Xenium | transcript point cloud + DAPI | YES (segmented) | SEGMENT + annotate -- do NOT deconvolve |
| MERFISH / MERSCOPE | subcellular | YES | SEGMENT + annotate -- do NOT deconvolve |
| CosMx SMI | subcellular | YES | SEGMENT + annotate -- do NOT deconvolve |

- DECONVOLVE: a capture array has no per-transcript cell assignment, so the spot is an unavoidable mixture and only mixture modeling recovers composition. The H&E underlay can count nuclei but cannot say which transcript came from which cell.
- SEGMENT (imaging platforms): the data are already single cells. The work is segmentation (assign transcripts to cells) then annotation (cluster + markers, or label-transfer). Running deconvolution here is conceptually WRONG -- it fabricates fractional mixtures inside cells that are already pure. See image-analysis for segmentation and single-cell/cell-annotation for label transfer.
- AMBIGUOUS (the near-single-cell middle): bins/beads still hold partial or multiple cells. When a high-quality registered image exists (Visium HD), morphology-driven cell RECONSTRUCTION (Bin2cell, StarDist/Cellpose nuclei expansion) is increasingly preferred over treating bins as fixed mixtures; without per-bead morphology (Slide-seqV2), assignment/deconvolution (often RCTD doublet-mode) remains standard. No settled consensus -- see high-resolution-binning.

## Governing Principle

A deconvolution result is a PROJECTION of the single-cell reference onto the spatial data. The reference is not a neutral input -- it is the dominant determinant of the answer, more than the algorithm.

The #1 trap is the missing or mismatched cell type. If a real type is ABSENT from the reference, its transcripts are forced onto the transcriptionally nearest type that IS present. The proportions still sum to 1 and look clean, and there is no internal warning -- garbage reference produces confident garbage proportions. The reference must match the TISSUE, the CONDITION (a healthy reference mis-estimates activated-immune or malignant states whose expression has shifted), and ideally the technology/protocol (snRNA-seq vs whole-cell dissociation bias propagates straight into the numbers). The reference is NOT ground truth.

Every benchmark reinforces the same lesson: reference quality and the cell-type abundance pattern swing results MORE than method choice, and a plain NNLS baseline outperforms nearly half the dedicated methods (Sang-Aram 2024 *eLife* 12:RP88431; Li 2022 *Nat Methods* 19:662-670). Method-shopping is the wrong lever; reference quality is the right one. Rare-type fractions (below a few percent) are the least reliable numbers in the output -- corroborate any rare type with its spatial marker genes before believing it. More tools is not more confidence: two NB-regression methods agreeing is pseudo-replication, not orthogonal validation. The defensible confidence move is reference-sensitivity analysis (perturb or swap the reference) plus an orthogonal modality.

Outputs are COMPOSITIONAL: per-spot proportions live on a simplex (sum to 1), so an increase in one type mechanically decreases the others. Downstream comparison must use CLR/ILR or compositional tests -- naive per-type t-tests/correlations on raw proportions are mis-specified.

## Choosing a method

The first axis is output: cell2location uniquely returns ABSOLUTE cell abundance (expected cell numbers per spot); the rest return proportions only. "Number of cells of type X" and "fraction of the spot that is type X" answer different biological questions. The second axis is runtime and whether spatial coordinates or a platform-shift correction are modeled. When competing methods exist, verify current best practice against the latest benchmark before committing -- this field moves fast.

| Method | Model | Reference | Output | Best when | Fails when |
|--------|-------|-----------|--------|-----------|------------|
| cell2location | Bayesian hierarchical NB (variational, pyro/scvi-tools) | yes | ABSOLUTE abundance | absolute counts wanted; large atlases; models platform shift | GPU-light setups (slow VI); over-trusting rare types |
| RCTD (spacexr) | Poisson + per-gene platform-effect random effect | yes | proportions | fast, widely used; doublet-mode for Slide-seq/high-res | non-R pipelines without rpy2 |
| stereoscope | Negative-binomial MLE of spot mixtures | yes | proportions | principled NB; in scvi-tools | speed (among slowest) |
| SPOTlight | Seeded NMF + NNLS | yes | proportions | fast, transparent, R/Bioconductor | mid-pack accuracy |
| SpatialDWLS | Dampened weighted least squares + marker enrichment | yes | proportions | fastest tier; consistently top in Li 2022 | needs Giotto |
| CARD | CAR-prior spatially-informed NMF regression | yes (can run ref-free) | proportions (smoothed) | spatially structured tissue; coordinates help | sharp composition boundaries (over-smooths) |
| Tangram | Deep-learning alignment (cell->voxel mapping) | yes | mapping (proportions as by-product) | transcript imputation; flexible platforms | pure proportion accuracy (it is a mapper) |
| STdeconvolve | Reference-FREE LDA topic model | NO | proportions + topic profiles | no matched reference exists; sanity check | types co-occurring in fixed ratios; topics need post-hoc annotation |

cell2location and RCTD are reliably top-tier across independent benchmarks (Li 2022; Sang-Aram 2024; *Nat Commun* 2023 14:1548). When no matched reference exists, STdeconvolve is the escape hatch -- it is immune to the missing-type trap because it uses no reference, but its topics are de novo and must be annotated afterward.

## cell2location step 1: reference signatures

**Goal:** Learn each cell type's per-gene expression signature from the annotated scRNA-seq reference, correcting for the technical/batch structure of the reference.

**Approach:** Filter genes to informative ones, fit a negative-binomial RegressionModel with the cell-type label and any batch as covariates, then export the posterior per-cluster mean expression. cell2location consumes RAW integer counts -- do NOT pass log-normalized data.

```python
import cell2location
import numpy as np
import scanpy as sc
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel

adata_ref = sc.read_h5ad('reference_scrna.h5ad')           # raw counts in .X
adata_ref.obs['cell_type'] = adata_ref.obs['cell_type'].astype('category')

selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
adata_ref = adata_ref[:, selected].copy()

# batch_key absorbs technical structure across reference samples; drop if a single batch
RegressionModel.setup_anndata(adata_ref, labels_key='cell_type', batch_key='sample')
mod = RegressionModel(adata_ref)
mod.train(max_epochs=250, accelerator='gpu')               # use_gpu= is deprecated; accelerator in {'gpu','cpu','auto'}
adata_ref = mod.export_posterior(adata_ref, sample_kwargs={'num_samples': 1000})

factors = adata_ref.uns['mod']['factor_names']
if 'means_per_cluster_mu_fg' in adata_ref.varm:
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{f}' for f in factors]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{f}' for f in factors]].copy()
inf_aver.columns = factors                                 # genes x cell_types signature matrix
```

## cell2location step 2: spatial mapping

**Goal:** Decompose each spatial spot into absolute cell-type abundances using the reference signatures.

**Approach:** Restrict both objects to shared genes, set up the Cell2location model with the signature matrix and the expected cells-per-spot, then train. N_cells_per_location is a tissue-dependent prior (Visium ~10-30, denser tissue higher); detection_alpha controls within-experiment normalization (20 is the tutorial default; raise toward 200 if technical variability in total counts is high).

```python
adata_vis = sc.read_h5ad('visium.h5ad')                    # raw counts
shared = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, shared].copy()
inf_aver = inf_aver.loc[shared, :]

cell2location.models.Cell2location.setup_anndata(adata_vis, batch_key='sample')
mod_sp = cell2location.models.Cell2location(adata_vis, cell_state_df=inf_aver,
                                            N_cells_per_location=30, detection_alpha=20)
mod_sp.train(max_epochs=30000, batch_size=None, train_size=1, accelerator='gpu')
adata_vis = mod_sp.export_posterior(adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod_sp.adata.n_obs})

# q05 = 5% posterior quantile = 'at least this many cells of this type are present' (conservative)
abundance = adata_vis.obsm['q05_cell_abundance_w_sf']      # ABSOLUTE expected cell numbers per spot
abundance.columns = factors
adata_vis.obs[abundance.columns] = abundance.values
```

cell2location returns absolute abundances. Convert to proportions only if proportions are the question -- doing so discards the information that distinguishes cell2location from the proportion-only methods.

## Handling compositional outputs

**Goal:** Compare cell-type composition across spots, regions, or conditions without the spurious negative correlations that the sum-to-1 constraint manufactures.

**Approach:** Map proportions out of the simplex with the centered log-ratio (CLR) before any correlation, t-test, or linear model. Add a small pseudocount because CLR is undefined at zero.

```python
import numpy as np

def clr(proportions, pseudocount=1e-6):
    p = proportions + pseudocount
    p = p / p.sum(axis=1, keepdims=True)
    log_p = np.log(p)
    return log_p - log_p.mean(axis=1, keepdims=True)       # subtract per-spot geometric-mean log

abund = adata_vis.obsm['q05_cell_abundance_w_sf'].values
proportions = abund / abund.sum(axis=1, keepdims=True)
clr_comp = clr(proportions)                                # now safe for downstream correlation / DE / t-tests
```

For differential abundance between conditions, prefer a compositional method (ALDEx2, scCODA, or a Dirichlet-multinomial model) over per-type t-tests on raw fractions; the same closed-data logic that breaks naive correlations breaks naive DA.

## Reference-free sanity check

**Goal:** Detect a missing-reference-type problem and recover composition when no matched scRNA-seq reference exists.

**Approach:** Run STdeconvolve (LDA topic model, R) on the spatial data alone, then check whether any de novo topic matches no reference cell type -- a topic with strong spatial structure and marker genes for a type absent from the reference is direct evidence the reference is incomplete.

```r
library(STdeconvolve)
counts <- cleanCounts(spatial_counts_matrix, min.lib.size = 100)   # genes x spots
corpus <- restrictCorpus(counts, removeAbove = 1.0, removeBelow = 0.05)
ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(8, 15))              # K = number of topics, scan a range
opt <- optimalModel(models = ldas, opt = 'min')
res <- getBetaTheta(opt, perc.filt = 0.05)
deconv_prop <- res$theta                                          # spots x topics proportions
# annotate topics post hoc via res$beta (topic gene profiles) against known markers
```

## Validating against marker genes

**Goal:** Confirm that an estimated cell-type fraction tracks the in-situ expression of that type's canonical markers, not the reference's wishful thinking.

**Approach:** For each type, correlate its estimated proportion across spots with the mean spatial expression of its marker genes; weak or negative correlation flags a misassignment or a reference mismatch.

```python
markers = {'T_cell': ['CD3D', 'CD3E', 'CD8A'], 'Macrophage': ['CD68', 'CD14', 'CSF1R'], 'Epithelial': ['EPCAM', 'KRT8']}
for ct, genes in markers.items():
    present = [g for g in genes if g in adata_vis.var_names]
    if not present or ct not in proportions_df.columns:
        continue
    expr = np.asarray(adata_vis[:, present].X.mean(axis=1)).ravel()
    corr = np.corrcoef(expr, proportions_df[ct].values)[0, 1]
    print(f'{ct}: marker-vs-proportion r = {corr:.3f}')        # low r -> suspect reference or misassignment
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Deconvolution "works" on Xenium/MERFISH/CosMx but fractions are nonsensical | Ran deconvolution on single-cell-resolution imaging data, inventing mixtures inside pure cells | Segment then annotate (image-analysis, single-cell/cell-annotation); deconvolution does not apply |
| A histologically present cell type appears in NO spot | Type is absent from the reference; its signal was silently reassigned to the nearest present type | Add the missing type to the reference; cross-check with STdeconvolve for an unmatched topic |
| Confident proportions, but disease/activated states look wrong | Condition mismatch -- healthy reference deconvolving diseased tissue | Use a condition-matched reference; activated/malignant states are transcriptionally far from healthy |
| Per-type t-test finds "everything changes oppositely" | Naive test on compositional (sum-to-1) proportions creates spurious negative correlation | CLR/ILR transform first; use ALDEx2/scCODA/Dirichlet-multinomial for differential abundance |
| Rare cell type fraction swings wildly between runs/references | Rare-type fractions are the least reliable output; abundance-pattern sensitivity | Treat <few-percent fractions skeptically; corroborate with spatial markers; do reference-sensitivity analysis |
| `ValueError`/garbage from cell2location after normalizing | Passed log-normalized data; cell2location needs RAW integer counts | Feed raw counts; stash normalized layers separately |
| `TypeError: unexpected keyword 'use_gpu'` | `use_gpu` deprecated in current scvi-tools | Use `accelerator='gpu'` (or 'cpu'/'auto') |
| Two methods agree, reported as validation | Two NB-regression methods are pseudo-replication, not orthogonal | Validate by perturbing the reference and against an orthogonal modality (matched imaging, in-situ markers) |
| Every spot contains a little of every immune type | Spot-edge transcript spillover / diffusion inflates apparent co-localization | RCTD doublet-mode or spillover-aware segmentation; treat ubiquitous low fractions with suspicion |

## Related Skills

- spatial-domains - group spots into tissue regions; a domain is a region, not a cell type or a niche
- high-resolution-binning - the AMBIGUOUS regime (Visium HD, Stereo-seq, Slide-seq): reconstruct cells vs deconvolve bins
- image-analysis - segment cells from imaging platforms, where deconvolution does not apply
- single-cell/cell-annotation - annotate the imaging cells after segmentation, and label the scRNA-seq reference
- single-cell/preprocessing - build a clean reference; the reference IS the result
- spatial-preprocessing - QC and normalize the spatial data before deconvolution
- spatial-visualization - map estimated proportions and abundances onto the tissue

## References

- Kleshchevnikov V, Shmatko A, Dann E, et al. (2022) Cell2location maps fine-grained cell types in spatial transcriptomics. Nature Biotechnology 40:661-671. DOI 10.1038/s41587-021-01139-4
- Cable DM, Murray E, Zou LS, et al. (2022) Robust decomposition of cell type mixtures in spatial transcriptomics (RCTD). Nature Biotechnology 40:517-526. DOI 10.1038/s41587-021-00830-w
- Andersson A, Bergenstrahle J, Asp M, et al. (2020) Single-cell and spatial transcriptomics enables probabilistic inference of cell type topography (stereoscope). Communications Biology 3:565. DOI 10.1038/s42003-020-01247-y
- Elosua-Bayes M, Nieto P, Mereu E, Gut I, Heyn H (2021) SPOTlight: seeded NMF regression to deconvolute spatial transcriptomics spots with single-cell transcriptomes. Nucleic Acids Research 49(9):e50. DOI 10.1093/nar/gkab043
- Biancalani T, Scalia G, Buffoni L, et al. (2021) Deep learning and alignment of spatially resolved single-cell transcriptomes with Tangram. Nature Methods 18(11):1352-1362. DOI 10.1038/s41592-021-01264-7
- Ma Y, Zhou X (2022) Spatially informed cell-type deconvolution for spatial transcriptomics (CARD). Nature Biotechnology 40:1349-1359. DOI 10.1038/s41587-022-01273-7
- Dong R, Yuan GC (2021) SpatialDWLS: accurate deconvolution of spatial transcriptomic data. Genome Biology 22:145. DOI 10.1186/s13059-021-02362-7
- Miller BF, Huang F, Atta L, Sahoo A, Fan J (2022) Reference-free cell type deconvolution of multi-cellular pixel-resolution spatially resolved transcriptomics data (STdeconvolve). Nature Communications 13:2339. DOI 10.1038/s41467-022-30033-z
- Sang-Aram C, Browaeys R, Seurinck R, Saeys Y (2024) Spotless, a reproducible pipeline for benchmarking cell type deconvolution in spatial transcriptomics. eLife 12:RP88431. DOI 10.7554/eLife.88431
- Li B, Zhang W, Guo C, et al. (2022) Benchmarking spatial and single-cell transcriptomics integration methods for transcript distribution prediction and cell type deconvolution. Nature Methods 19:662-670. DOI 10.1038/s41592-022-01480-9
