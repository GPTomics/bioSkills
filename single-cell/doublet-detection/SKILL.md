---
name: bio-single-cell-doublet-detection
description: Detect and remove doublets (two or more cells in one droplet) from single-cell RNA-seq using scDblFinder (R), Scrublet (Python), and DoubletFinder (R). Use when flagging artificial intermediate populations before clustering, setting the expected doublet rate from recovered-cell counts, running detection per sample before integration, choosing between simulate-and-score methods, or interpreting a non-bimodal score histogram.
tool_type: mixed
primary_tool: scDblFinder
---

## Version Compatibility

Reference examples tested with: scanpy 1.10+, scDblFinder 1.16+, Seurat 5.0+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Doublet Detection

**"Remove doublets from my data"** -> Flag droplets that captured two or more cells, which masquerade as fake intermediate cell states.
- Python: `sc.pp.scrublet()` per sample on raw counts
- R: `scDblFinder(sce, samples=...)` per sample on raw counts

## Governing Principle

Doublets fabricate fake biology, so the goal is not a "doublet-free" dataset but avoiding false conclusions. Three facts govern every decision.

Doublets create fake intermediate populations. A heterotypic doublet (two distinct types, e.g. T cell + monocyte) sums to a profile that lands between clusters and reads as a novel "transitional" state - the most damaging failure mode, because it corrupts trajectory inference and RNA velocity by building false bridges between lineages. Treat any small cluster co-expressing two lineage programs (CD3+LYZ, EPCAM+PTPRC) as doublet-suspect until proven otherwise.

Detect per sample, before integration or clustering. A doublet is a physical event within one droplet in one capture, so two cells from different samples can never share one - any cross-sample doublet called on a merged object is meaningless. Merging also corrupts the kNN/PCA neighborhood that scoring depends on. scDblFinder's `samples=` handles this internally; Scrublet and DoubletFinder must be looped per sample. All three want raw counts after basic QC.

Removal is never complete, and over-removal deletes real cells. Homotypic doublets (two cells of the same type) sum to a profile that looks like one bigger cell of that type and are nearly invisible to any expression-based method, so reported "doublet rates" only cover the heterotypic-detectable fraction. Conversely, doublet scores correlate with total counts, the same axis as count-based QC, so aggressive filtering on both double-penalizes and strips genuine high-RNA populations (megakaryocytes, plasma cells, large neurons). Coordinate the two filters and prefer flag-and-inspect over blind deletion.

## Expected Doublet Rate

10X Chromium loading is near-Poisson, so the multiplet rate scales roughly linearly with recovered cells: ~0.8% per 1,000 cells recovered (`dbr.per1k = 0.008`).

| Cells recovered (~) | Expected rate (~) |
|---------------------|-------------------|
| 1,000 | 0.8% |
| 2,000 | 1.6% |
| 5,000 | 3.9% |
| 10,000 | 7.6-8% |

Rule: rate ~= 0.008 x recovered/1000. Always set the expected rate from the actual recovered-cell count of that lane; Scrublet's flat `expected_doublet_rate=0.05` is a placeholder, not a recommendation. High-throughput chips have lower per-cell rates. For multiplexed pools (genotype/HTO-demultiplexed), the physical doublet rate is set by TOTAL lane loading, not the demultiplexed subset: deriving the rate from one sample's cells underestimates it (four 5k samples in one 20k lane is ~15% real, not the ~3.9% implied by 5k), so set the rate from the total lane cell count.

## Heterotypic, Homotypic, Neotypic

| Type | Composition | Detectability |
|------|-------------|---------------|
| Heterotypic | Two transcriptionally distinct types | Detectable; lands between clusters; the dangerous "fake transitional" ones |
| Homotypic | Two cells of the same type | Nearly undetectable by expression; persists after removal |
| Neotypic | Heterotypic blend occupying a region no singlet occupies | Most detectable; most misleading if missed (looks like a rare new type) |

`modelHomotypic`-style adjustments only change the number expected to be detectable; they cannot recover undetectable homotypic doublets.

## Choosing a Method

| Method | Model | Use when | Fails / weak when |
|--------|-------|----------|-------------------|
| scDblFinder (R) | xgboost on kNN features vs simulated doublets | Default; best accuracy-speed balance; built-in per-sample via `samples=` | R/Bioconductor only |
| Scrublet (Python) | kNN density of simulated doublets | scanpy-native pipelines | Auto-threshold fails on unimodal histograms; loop per sample manually |
| DoubletFinder (R) | pANN from PC neighborhood | Legacy Seurat workflows | Brittle; `pK` needs per-dataset sweep; `*_v3` names removed; Seurat-version-coupled |
| solo (Python) | scVI VAE + classifier | Have a trained scVI model; GPU available | Heavier setup |
| scds cxds/bcds/hybrid (R) | Co-expression / boosted tree | Fast first pass on very large data | Lower accuracy than scDblFinder/DoubletFinder |

scDblFinder is the 2024-2026 best-balance default and is recommended by sc-best-practices. The older Xi and Li 2021 ranking ("DoubletFinder is most accurate") predates major scDblFinder improvements and is superseded - do not cite it against current scDblFinder. Methods compete and drift; verify current standing against the installed tool's docs before committing.

## scDblFinder (R, recommended)

**Goal:** Call doublets with a fast gradient-boosted classifier, per sample, with the rate inferred from cell count.

**Approach:** Convert to SingleCellExperiment, pass the per-sample key so each capture is processed independently, then read the class/score back.

```r
library(scDblFinder)
library(SingleCellExperiment)

sce <- as.SingleCellExperiment(seurat_obj)                 # or build directly from a counts matrix
sce <- scDblFinder(sce, samples = 'sample_id')             # per-capture; dbr defaults from cell count via dbr.per1k=0.008
table(sce$scDblFinder.class)                               # adds scDblFinder.class ('singlet'/'doublet') and .score
seurat_obj$scDblFinder_class <- sce$scDblFinder.class
seurat_obj$scDblFinder_score <- sce$scDblFinder.score
```

`clusters=NULL` (default) generates purely random artificial doublets and is generally recommended; pass a vector for cluster-based generation. `dbr=NULL` computes the rate from cell count; set `dbr`/`dbr.sd` explicitly to encode a known loading.

## Scrublet (Python)

**Goal:** Score doublets in a scanpy pipeline, per sample, with the rate set from recovered cells.

**Approach:** Use the maintained `sc.pp.scrublet` path on raw counts; set `expected_doublet_rate` per lane; inspect the histogram when the auto-threshold looks wrong.

```python
import scanpy as sc

n_cells = adata.n_obs
expected_rate = 0.008 * n_cells / 1000                     # from recovered cells, not the 0.05 placeholder
sc.pp.scrublet(adata, expected_doublet_rate=expected_rate)  # adds obs['doublet_score'], obs['predicted_doublet']
# auto-threshold needs a bimodal histogram; if unimodal, inspect uns['scrublet'] and set threshold manually
adata_singlets = adata[~adata.obs['predicted_doublet']].copy()
```

For pooled samples, loop `sc.pp.scrublet(adata[adata.obs.sample == s], ...)` per sample (or pass `batch_key`), never on the merged object.

## DoubletFinder (R, legacy Seurat)

**Goal:** Run DoubletFinder on a fully preprocessed Seurat object, tuning `pK` and homotypic-adjusting the expected count.

**Approach:** Sweep `pK`, pick the BCmvn maximum, then set `nExp` from the rate adjusted for the homotypic fraction.

```r
library(DoubletFinder)                                     # *_v3 function names were removed in Nov 2023; verify installed API

sweep.res <- paramSweep(seurat_obj, PCs = 1:20, sct = FALSE)
bcmvn <- find.pK(summarizeSweep(sweep.res, GT = FALSE))
pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))   # no default pK; tune per dataset

rate <- 0.008 * ncol(seurat_obj) / 1000
nExp <- round(rate * ncol(seurat_obj))
nExp <- round(nExp * (1 - modelHomotypic(seurat_obj$seurat_clusters)))  # discount undetectable homotypic doublets
seurat_obj <- doubletFinder(seurat_obj, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp, sct = FALSE)
```

`pN` (artificial-doublet proportion) defaults to 0.25 and performance is largely insensitive to it. DoubletFinder requires a normalized, PCA'd, clustered object and is the most version-sensitive of the three.

## Deeper Cautions

Simulated doublets are a model, not the real thing: real doublets share one RT/PCR reaction (capture competition, barcode effects), so simulated-doublet density only approximates where real doublets sit, and even the best method has a low ceiling (max mean AUPRC ~0.537 in Xi and Li 2021 - every method misses a lot). Over-removal culls proliferating (S/G2M) and genuine transitional cells that legitimately score high, so cross-check removed cells against cell-cycle and activation signatures. When available, experimental ground truth beats inference: cell hashing (CITE-seq HTOs) and MULTI-seq call inter-sample doublets directly regardless of expression similarity (catching even cross-sample homotypic doublets), and serve as a complementary filter. Heavy ambient RNA can mimic co-expression and nudge scores, so handle empty droplets and ambient RNA first (see single-cell/preprocessing).

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Doublet calls look random / too many | Run on merged multi-sample data | Run per sample before integration (`samples=` or loop) |
| Auto-threshold splits the histogram badly | Scrublet histogram is unimodal | Inspect the histogram and set `threshold` manually |
| A high-RNA cell type was wiped out | Count-based QC and doublet removal double-penalized the same axis | Coordinate the filters; do not stack aggressive cutoffs |
| "Novel transitional state" co-expresses two lineages | Heterotypic doublets masquerading as a cluster | Confirm per-sample detection; check marker co-expression / hashing before claiming a new type |
| Trajectory has an implausible bridge between lineages | Doublets forming a false intermediate | Remove/flag doublets before trajectory inference |
| Reported "0% doublets" | Homotypic doublets are invisible | Do not claim doublet-free; report only the detectable fraction |
| DoubletFinder call errors after a Seurat upgrade | `*_v3` names removed; API drift | Use current function names; re-tune `pK` |
| Expected rate clearly wrong | Used a package default | Set rate from recovered cells (~0.008 x cells/1000) |
| Multiplexed pool underestimates doublets | Rate derived from one demultiplexed sample, not total lane | Set the expected rate from total capture-lane cells |

## Related Skills

- single-cell/preprocessing - QC and ambient-RNA handling before doublet detection
- single-cell/hashing-demultiplexing - Hashtag-based cross-sample doublet calling (complements expression-based detection)
- single-cell/data-io - load raw per-sample matrices before processing
- single-cell/clustering - run clustering after doublet removal
- single-cell/batch-integration - integrate samples only after per-sample doublet calling
- single-cell/trajectory-inference - doublets create false bridges; remove them first

## References

- Wolock SL, Lopez R, Klein AM (2019) Scrublet: computational identification of cell doublets in single-cell transcriptomic data. Cell Systems 8(4):281-291.e9. DOI 10.1016/j.cels.2018.11.005
- McGinnis CS, Murrow LM, Gartner ZJ (2019) DoubletFinder: doublet detection in single-cell RNA sequencing data using artificial nearest neighbors. Cell Systems 8(4):329-337.e4. DOI 10.1016/j.cels.2019.03.003
- Germain P-L, Lun A, Macnair W, Robinson MD (2021) Doublet identification in single-cell sequencing data using scDblFinder. F1000Research 10:979. DOI 10.12688/f1000research.73600
- Xi NM, Li JJ (2021) Benchmarking computational doublet-detection methods for single-cell RNA sequencing data. Cell Systems 12(2):176-194.e6. DOI 10.1016/j.cels.2020.11.008
- Bernstein NJ, Fong NL, Lam I, et al. (2020) Solo: doublet identification in single-cell RNA-seq via semi-supervised deep learning. Cell Systems 11(1):95-101.e5. DOI 10.1016/j.cels.2020.05.010
- Bais AS, Kostka D (2020) scds: computational annotation of doublets in single-cell RNA sequencing data. Bioinformatics 36(4):1150-1158. DOI 10.1093/bioinformatics/btz698
- McGinnis CS, Patterson DM, Winkler J, et al. (2019) MULTI-seq: sample multiplexing for single-cell RNA sequencing using lipid-tagged indices. Nature Methods 16(7):619-626. DOI 10.1038/s41592-019-0433-8
- Heumos L, Schaar AC, Lance C, et al. (2023) Best practices for single-cell analysis across modalities. Nature Reviews Genetics 24:550-572. DOI 10.1038/s41576-023-00586-w
