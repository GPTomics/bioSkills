---
name: bio-spatial-transcriptomics-spatial-preprocessing
description: Quality control, filtering, and normalization for spatial transcriptomics (Visium, Visium HD, Xenium, MERFISH/MERSCOPE, CosMx, Slide-seq) with Squidpy and Scanpy. Use when setting QC floors that do NOT delete real low-count imaging cells (an scRNA min_counts=500 floor deletes nearly every Xenium cell, whose vector is tens-to-low-hundreds of transcripts); deciding whether to normalize at all when library size carries spatial biology rather than pure technical depth; choosing cell-volume/area normalization over Pearson residuals for skewed targeted panels; reading negative-control-probe / blank-barcode false-discovery rates; and inspecting QC spatially on the tissue rather than only in violins.
tool_type: python
primary_tool: squidpy
---

## Version Compatibility

Reference examples tested with: squidpy 1.5+, scanpy 1.10+, anndata 0.10+, spatialdata 0.2+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Spatial Preprocessing

**"QC and normalize my spatial data"** -> Flag and remove low-quality spots/cells, then put counts on a scale fit for downstream domain and marker analysis -- but the right QC floors and the right normalization both depend on which side of the platform fork the data sits.
- Sequencing/spot (Visium, Visium HD, Slide-seq, Stereo-seq): a spot/bin is mini-bulk over 1-10 cells; QC on UMI/spot, genes/spot, mito-%, cells/spot; normalization must respect that library size tracks cellularity.
- Imaging/in-situ (Xenium, MERSCOPE/MERFISH, CosMx, seqFISH): a cell is segmentation-derived, carries tens-to-low-hundreds of transcripts over a TARGETED panel; QC on a low transcript floor, cell area, and negative-control FDR; normalization must not be gene-count-based.

## The platform-class fork (decide this first)

The first question on any spatial dataset is which assay family produced it, because it changes every QC threshold and the entire normalization decision.

| Axis | Sequencing/spot (Visium, Slide-seq, Stereo-seq) | Imaging/in-situ (Xenium, MERSCOPE, CosMx) |
|---|---|---|
| Unit | spot/bin = 1-10-cell mixture | segmentation-derived single cell |
| Counts/unit | hundreds-thousands UMI | tens-low hundreds transcripts |
| Gene space | whole-transcriptome (poly-A) or probe panel | TARGETED panel (100-1000), genes/cell ceilinged at panel size |
| Mito-% QC | available (mixed-cell average) | usually impossible (mito off-panel) |
| Specificity metric | none native | negative-control-probe / blank-barcode FDR |
| Library-size meaning | confounds cells-per-spot + cellularity | confounds cell SIZE/AREA + segmentation error |

## Governing Principle

In single-cell RNA-seq library size is a technical nuisance to divide out. In spatial transcriptomics LIBRARY SIZE CARRIES BIOLOGY, and that single fact governs both QC and normalization. On Visium, total counts per spot are spatially structured and correlated with anatomy because they confound with the number of cells per spot and tissue cellularity (Bhuva 2024 *Genome Biol* 25:99). On imaging platforms, total counts per cell confound with cell SIZE/AREA -- a physically larger segmented cell holds more molecules for purely geometric reasons -- and with segmentation error itself. Naively dividing library size out (CP10k, log1p, scran pooling) therefore removes real spatially-structured biology and measurably degrades spatial-domain detection.

For imaging the bias is UPSTREAM of any residual model. Because a targeted panel is small, hand-curated, and skewed toward a few high markers, any gene-count-based size factor is dominated by a handful of genes and becomes panel-composition-dependent. Atta and Fan 2024 (*Genome Biol* 25:153) compared library-size, Pearson/SCTransform, DESeq2, TMM, and volume/area normalization on skewed panels: the four gene-count methods inject region-specific bias of up to ~13% DE error and fold-change SIGN REVERSAL in up to 19% of genes, while volume/area normalization avoids it because its denominator is independent of panel composition. The load-bearing consequence: Pearson residuals and SCTransform, the gold standard for whole-transcriptome scRNA, do NOT rescue a skewed imaging panel -- they are still gene-count-based and the bias sits in the panel design. The fix is non-gene-count normalization (cell volume/area, Moffitt-style) or spatially-aware joint modeling (SpaNorm).

A clean violin plot proves nothing. A count floor that silently deletes a spatial cluster of small cells (lymphocytes), a normalization that erases a cellularity gradient, and a focal hybridization failure all survive a tidy genes-per-cell distribution. The dangerous artifacts are spatial, so QC must be inspected spatially.

## scRNA QC thresholds are wrong for imaging

Carrying scRNA defaults onto imaging data deletes the data. Imaging cells carry tens-to-low-hundreds of transcripts -- one to two orders of magnitude below droplet scRNA -- so an `min_counts=500` floor removes nearly every real cell. Genes/cell can never exceed the panel size, so the "high genes = doublet" heuristic is meaningless (the ceiling is the panel, not a doublet). Mito genes are usually off-panel, so `pct_counts_mt` QC is often impossible. Worst, an aggressive count floor preferentially deletes the smallest REAL cells (lymphocytes, neutrophils), biasing tissue composition rather than removing noise.

### QC-metric-by-platform table

| Metric | Visium (spot) | Imaging (Xenium/MERSCOPE/CosMx) | Rationale / trap |
|---|---|---|---|
| counts/unit | UMI/spot; no universal floor; OSTA DLPFC flags <600 | transcripts/cell tens-low hundreds; community floor ~10 (Squidpy), some 20 | scRNA `min_counts=500` deletes nearly every imaging cell; floor confounds cellularity/cell-size |
| genes/unit | genes/spot; OSTA flags <400 | genes/cell CEILING = panel size (100-1000) | "high genes = doublet" meaningless for imaging |
| mito-% | mixed-cell average; OSTA flags >0.28 (brain) | usually off-panel -> impossible | tissue-dependent; brain tolerates higher |
| cell area (um^2) | n/a (spot is fixed) | MAD-based on counts/area | flags over/under-segmentation; no fixed vendor min |
| negative-control FDR | n/a | THE imaging specificity metric | false-discovery proxy; no scRNA analogue |
| cells/spot | nuclei estimate; OSTA flags >10 | n/a | confirms spot is a mixture |

Thresholds are tissue-dependent and "somewhat arbitrary" (the OSTA Visium worked example flags UMI<600, genes<400, mito>0.28, cells/spot>10 on DLPFC, removing 32/3639 spots -- a starting point, not a law). The imaging floor of ~10 transcripts/cell is a Squidpy/community convention, NOT a vendor specification; community CosMx floors run higher (commonly ~20 counts/cell, scaling up with plex), so confirm the cutoff against the panel and tissue rather than copying a number.

**Goal:** Annotate negative controls and mito genes (where present), compute QC metrics, and set platform-appropriate floors without deleting real low-count cells.

**Approach:** Branch on the fork. For imaging, identify control-probe prefixes, filter on a low transcript floor and cell area; for spot data, use UMI/genes/mito floors. Always compute metrics, then look at them spatially before cutting.

```python
import squidpy as sq
import scanpy as sc
import numpy as np

# Imaging branch: control features carry platform-specific prefixes -- they are the specificity ruler, not genes
ctrl_prefixes = ('NegControlProbe', 'NegControlCodeword', 'BLANK', 'Blank', 'NegPrb')   # Xenium / MERFISH / CosMx
adata.var['control'] = adata.var_names.str.startswith(ctrl_prefixes)
adata.var['mt'] = adata.var_names.str.startswith(('MT-', 'mt-'))                          # usually empty on imaging panels
sc.pp.calculate_qc_metrics(adata, qc_vars=['control', 'mt'], percent_top=None, inplace=True)
# inplace defaults to False and returns DataFrames; pass inplace=True to write .obs/.var
```

## Negative controls -- the imaging specificity metric

Imaging platforms include features that decode to nothing biological: Xenium negative-control PROBES (off-target binding) plus negative-control CODEWORDS (pure optical/decoding error), MERFISH/MERSCOPE blank barcodes (valid codewords with no probe), CosMx NegPrb (alien synthetic sequences). They are the only native false-discovery proxy in spatial data. The canonical metric is FDR = mean counts per control feature / mean counts per real gene; the community-acceptable band is roughly <=1-5% of signal (Xenium typically <0.1%, MERFISH ~4%, CosMx highest). Compute it before trusting any gene-level claim, and treat controls as a panel-wide QC gate -- not as genes to cluster on.

**Goal:** Quantify the per-feature false-discovery rate and drop controls before normalization and clustering.

**Approach:** Average per-feature counts within the control set and within real genes, take the ratio, then subset the matrix to real genes only.

```python
ctrl = adata.var['control'].values
mean_ctrl = np.asarray(adata[:, ctrl].X.sum(axis=0)).ravel().mean() if ctrl.any() else 0.0
mean_gene = np.asarray(adata[:, ~ctrl].X.sum(axis=0)).ravel().mean()
fdr = mean_ctrl / mean_gene if mean_gene else float('nan')
print(f'negative-control FDR: {fdr:.4f}  (band ~<=0.01-0.05)')
adata = adata[:, ~ctrl].copy()   # controls are a QC ruler, never clustering features
```

## Inspect QC spatially, then filter

A QC gradient across the section -- counts falling toward one edge, mito rising in a corner -- is a technical artifact (edge effects, uneven permeabilization, focal hybridization failure), not biology, and a violin plot hides it. Always map QC onto tissue coordinates before choosing thresholds, and check WHERE the cells slated for removal actually fall.

**Goal:** Reveal spatially-structured quality artifacts and confirm a proposed floor is not removing a coherent tissue region.

**Approach:** Color the spatial scatter by each QC metric; a smooth spatial gradient signals a technical artifact to address (or model) rather than threshold away.

```python
sq.pl.spatial_scatter(adata, color=['total_counts', 'n_genes_by_counts'], shape=None, ncols=2)
# shape=None renders points (imaging/Slide-seq); omit it for Visium hex spots with a tissue image
```

**Goal:** Apply platform-appropriate floors that remove debris and segmentation failures without biasing composition.

**Approach:** Imaging -- low transcript floor plus a cell-area sanity bound. Spot -- UMI/genes/mito floors. Either way, filter genes seen in too few units last.

```python
# Imaging floor: ~10 transcripts/cell is a Squidpy/community convention, NOT a vendor spec
sc.pp.filter_cells(adata, min_counts=10)
if 'cell_area' in adata.obs:
    lo, hi = adata.obs['cell_area'].quantile([0.01, 0.99])     # trim segmentation over/under-calls, tissue-dependent
    adata = adata[(adata.obs['cell_area'] > lo) & (adata.obs['cell_area'] < hi)].copy()
sc.pp.filter_genes(adata, min_cells=5)

# Spot branch instead (Visium): tissue-dependent floors -- the OSTA DLPFC example, not universal law
# sc.pp.filter_cells(adata, min_counts=600)
# sc.pp.filter_cells(adata, min_genes=400)
# adata = adata[adata.obs['pct_counts_mt'] < 28].copy()
```

## Normalization -- the central decision

Do not reach reflexively for `normalize_total` + `log1p`. The shipped Squidpy tutorials run it for both Xenium and MERFISH, so it is the de-facto default -- and it is exactly what the benchmark papers argue is biased for spatial data. Decide deliberately from the table, and because methods compete here, verify current best practice against the installed tool's docs and the latest benchmarks before committing.

### Normalization-method table

| Method | Assumption | Best when | Fails when |
|---|---|---|---|
| `normalize_total` + `log1p` | library size = pure technical nuisance | cross-platform comparability; quick default; tool tutorials | spatial -- removes spatially-structured biology (Bhuva 2024); imaging skewed panel |
| Analytic Pearson residuals | closed-form NB offset; depth as fixed offset | whole-transcriptome Visium HVG/PCA | imaging targeted panel -- still gene-count-based, inherits panel-skew bias (Atta/Fan 2024) |
| SCTransform v2 | regularized NB GLM, depth slope fixed | whole-transcriptome UMI / Visium | does NOT fix skewed imaging panels (gene-count-based) |
| Cell volume/area | concentration is the biological quantity | IMAGING skewed panel (Moffitt-style) | denominator needs reliable segmentation; cannot fix segmentation error itself |
| SpaNorm (spatially-aware) | library size and biology are entangled; remove only library-size component | spot AND imaging; preserve spatial structure | newer; R/Bioconductor |

**Goal (spot, whole-transcriptome):** Stabilize depth for HVG/PCA while keeping raw counts, accepting that crude library-size division can blur domains.

**Approach:** Stash raw counts, then either run the standard log1p pipeline knowingly or prefer analytic Pearson residuals for feature selection on whole-transcriptome Visium.

```python
adata.layers['counts'] = adata.X.copy()                    # stash raw -- HVG flavors and re-normalization need it
sc.pp.normalize_total(adata)                               # target_sum=None scales to the dataset MEDIAN, not the arbitrary 1e4
sc.pp.log1p(adata)                                         # library size carries biology -- this can blur spatial domains
# Whole-transcriptome Visium feature selection alternative (gene-count-based, fine here, NOT for imaging panels):
# sc.experimental.pp.normalize_pearson_residuals(adata)
```

**Goal (imaging, targeted panel):** Normalize without injecting panel-composition bias, using a denominator independent of gene counts.

**Approach:** Divide each cell's counts by its segmented area/volume (copies per unit area), then log-transform -- the Moffitt-style fix benchmarks favour over Pearson residuals for skewed panels.

```python
adata.layers['counts'] = adata.X.copy()
if 'cell_area' in adata.obs:                               # area/volume denominator is panel-composition-independent
    sf = adata.obs['cell_area'].values / adata.obs['cell_area'].median()
    adata.X = adata.X / sf[:, None]
    sc.pp.log1p(adata)
# If no segmentation area is available, prefer SpaNorm (R) over reflexive normalize_total on a skewed panel
```

## Common Errors

| Symptom | Cause | Fix |
|---|---|---|
| Nearly all imaging cells filtered out | scRNA `min_counts=500` floor on tens-of-transcript cells | Use a low floor (~10 transcripts/cell); branch QC on the platform fork |
| Spatial domains blur / merge after normalization | Crude library-size division erased a cellularity/anatomy gradient (library size carries biology) | Prefer SpaNorm or volume/area; if using log1p, know it can blur domains |
| Fold-change sign flips between normalizations | Gene-count size factor is panel-composition-dependent on a skewed imaging panel | Use non-gene-count (cell area/volume) normalization; Pearson/SCT do NOT fix it |
| `pct_counts_mt` is all zero / NaN | Mito genes are not on the targeted imaging panel | Skip mito-% QC for imaging; QC on transcript floor + cell area instead |
| "Doublet" cells flagged by high gene count | genes/cell ceiling IS the panel size -- not a doublet signal | Drop the high-genes heuristic for imaging; use cell area / spatial doublets |
| Control features cluster as their own group | Negative-control probes/codewords left in the matrix | Compute control FDR, then subset to real genes before clustering |
| Smallest cell type vanished after filtering | A count floor preferentially deleted small real cells (lymphocytes) | Inspect spatially where cuts fall; lower the floor; check composition before/after |
| Claimed a novel cell-state signature from imaging marker genes | An imaging panel (even 5,000-plex) is pre-selected for KNOWN biology; off-panel genes are absent by design, not by expression | Treat absence of an off-panel gene as uninformative; de-novo state discovery is bounded by the panel -- corroborate on whole-transcriptome data before claiming novelty |
| QC looks fine in violins but a region is empty | Spatial QC gradient (edge/permeabilization artifact) invisible in violins | Map QC onto tissue with `sq.pl.spatial_scatter` before thresholding |
| Counts inflated ~2x after re-running normalization | Normalized already-normalized data | Normalize raw once; restore from `layers['counts']` |

## Related Skills

- spatial-data-io - load Visium/Xenium/MERFISH and reach the molecule table vs the segmentation-derived matrix
- image-analysis - segment cells from imaging data, the upstream error source that sets imaging QC and cell area
- spatial-deconvolution - the next step for spot data, where library size and reference choice decide proportions
- single-cell/preprocessing - the scRNA QC/normalization baseline these thresholds deliberately depart from
- single-cell/clustering - cluster the QC'd cells; resolution is not a truth knob
- single-cell/cell-annotation - label-transfer typing for a targeted panel (de-novo marker discovery is panel-bounded)

## References

- Bhuva DD, Tan CW, Salim A, et al. (2024) Library size confounds biology in spatial transcriptomics data. Genome Biology 25:99. DOI 10.1186/s13059-024-03241-7
- Atta L, Clifton K, Anant M, Aihara G, Fan J (2024) Gene count normalization in single-cell imaging-based spatially resolved transcriptomics. Genome Biology 25:153. DOI 10.1186/s13059-024-03303-w
- Lause J, Berens P, Kobak D (2021) Analytic Pearson residuals for normalization of single-cell RNA-seq UMI data. Genome Biology 22:258. DOI 10.1186/s13059-021-02451-7
- Palla G, Spitzer H, Klein M, et al. (2022) Squidpy: a scalable framework for spatial omics analysis. Nature Methods 19:171-178. DOI 10.1038/s41592-021-01358-2
- Salim A, Bhuva DD, Chen C, et al. (2025) SpaNorm: spatially-aware normalisation for spatial transcriptomics data. Genome Biology 26:109. DOI 10.1186/s13059-025-03565-y
- Moffitt JR, Bambah-Mukku D, Eichhorn SW, et al. (2018) Molecular, spatial, and functional single-cell profiling of the hypothalamic preoptic region. Science 362:eaau5324. DOI 10.1126/science.aau5324
- Janesick A, Shelansky R, Gottscho AD, et al. (2023) High resolution mapping of the tumor microenvironment using integrated single-cell, spatial and in situ analysis. Nature Communications 14:8353. DOI 10.1038/s41467-023-43458-x
- Maynard KR, Collado-Torres L, Weber LM, et al. (2021) Transcriptome-scale spatial gene expression in the human dorsolateral prefrontal cortex. Nature Neuroscience 24:425-436. DOI 10.1038/s41593-020-00787-0
