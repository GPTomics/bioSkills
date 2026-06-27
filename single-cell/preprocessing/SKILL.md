---
name: bio-single-cell-preprocessing
description: Quality control, ambient-RNA handling, normalization, and feature selection for single-cell RNA-seq using Scanpy (Python) and Seurat (R). Use when filtering low-quality cells with MAD-adaptive thresholds, setting tissue-aware mito cutoffs, removing ambient RNA (SoupX/CellBender/DecontX), choosing a normalization (shifted-log vs scran vs sctransform vs Pearson residuals), selecting highly variable genes, or deciding whether to scale and regress out covariates.
tool_type: mixed
primary_tool: Seurat
---

## Version Compatibility

Reference examples tested with: scanpy 1.10+, Seurat 5.0+, scran 1.30+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Single-Cell Preprocessing

**"Preprocess my scRNA-seq data"** -> Remove bad barcodes, correct technical biases, and select informative genes before dimensionality reduction.
- Python: `calculate_qc_metrics()` -> filter -> `normalize_total()`+`log1p()` -> `highly_variable_genes()`
- R: QC -> `NormalizeData()` or `SCTransform()` -> `FindVariableFeatures()` -> `ScaleData()`

## Governing Principle

Every preprocessing choice propagates to every downstream result, and the two highest-leverage choices both encode a hidden biological assumption.

Normalization assumes near-constant total mRNA per cell. Shifted-log CP10k, scran deconvolution, and sctransform all divide by a per-cell size factor meant to capture only capture-efficiency and sequencing depth. They silently assume total transcriptome size is roughly constant across cell types, so a cell's total UMI count is a pure technical nuisance. This is false for plasma/antibody-secreting cells, secretory epithelia, hepatocytes, large neurons, and S/G2M cells, which carry 2-10x more mRNA. Dividing them to a common total deflates every gene that is not one of their few dominant transcripts (a compositional see-saw), and partially erases real proliferation biology. The honest framing: single-cell measures proportions, not absolute amounts.

QC metrics are biology metrics in disguise. `pct_counts_mt` conflates apoptosis, genuine metabolic activity (cardiomyocytes/hepatocytes/muscle are constitutively 20-40% mito and healthy), dissociation stress, and technical contamination. A flat global mito cutoff deletes entire healthy parenchymal populations and the survivors still cluster cleanly, so the loss is invisible. Use adaptive, tissue-aware thresholds and treat all three QC covariates jointly.

A beautiful UMAP proves nothing. Compositional normalization bias, deleted high-mito parenchyma, dissociation-stress clusters, ambient-induced co-expression, and residual homotypic doublets are all compatible with tidy clusters. The dangerous artifacts are precisely the ones that do not look like artifacts.

## Canonical Pipeline Order

1. Load the RAW (unfiltered) droplet matrix.
2. Empty-droplet calling (EmptyDrops, FDR<=0.001 on raw) or CellBender (folds calling + denoising).
3. Ambient-RNA removal (optional; SoupX/DecontX/CellBender) - BEFORE QC, because it needs the soup estimate.
4. QC filtering: cells (MAD on counts/genes/mito) + genes (`min_cells`).
5. Doublet detection - per sample, on raw counts (see single-cell/doublet-detection).
6. Normalization - shifted-log default, or scran/Pearson.
7. HVG selection - mind the raw-vs-lognorm input per flavor.
8. Scaling (optional, increasingly skipped).
9. PCA on HVG (~50 comps), then neighbors/clustering.

Ambient correction needs the raw matrix and must precede QC filtering; once subset to cells, the soup estimate is gone. Doublet detection runs on raw counts, so stash counts before normalizing.

Steps 1-5 are per-sample operations performed BEFORE merge or integration: empty-droplet calling, ambient removal (SoupX `load10X` is inherently per-run), adaptive QC, and doublet detection all reason about one capture's droplet population, so QC-then-merge is correct and merge-then-QC leaks batch effects into every threshold and contaminates the soup and doublet-scoring neighborhoods. Merge only after each sample is cleaned.

## Quality Control

**Goal:** Remove empty/dying/stressed barcodes using data-driven thresholds that do not delete real cell types.

**Approach:** Annotate mito/ribo/hemoglobin gene sets, compute joint QC metrics, then flag outliers by median absolute deviation (MAD) on the log scale rather than fixed cutoffs.

```python
import scanpy as sc
import numpy as np
from scipy.stats import median_abs_deviation

adata.var['mt'] = adata.var_names.str.startswith('MT-')                       # mouse: 'mt-'
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
adata.var['hb'] = adata.var_names.str.contains(r'^HB[ABDEGMQZ]\d*(?!\w)')      # explicit subunits, not legacy ^HB[^(P)]
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo', 'hb'], percent_top=[20], log1p=True, inplace=True)
# inplace defaults to False and returns DataFrames; pass inplace=True to write .obs/.var

def is_outlier(adata, metric, nmads):
    M = adata.obs[metric]
    return (M < np.median(M) - nmads * median_abs_deviation(M)) | (np.median(M) + nmads * median_abs_deviation(M) < M)

adata.obs['outlier'] = (is_outlier(adata, 'log1p_total_counts', 5) | is_outlier(adata, 'log1p_n_genes_by_counts', 5)
                        | is_outlier(adata, 'pct_counts_in_top_20_genes', 5))
adata.obs['mt_outlier'] = is_outlier(adata, 'pct_counts_mt', 3) | (adata.obs['pct_counts_mt'] > 8)
adata = adata[~(adata.obs['outlier'] | adata.obs['mt_outlier'])].copy()
sc.pp.filter_genes(adata, min_cells=3)
```

When samples differ in depth/quality or were sequenced in separate batches, compute MAD thresholds PER SAMPLE, not globally: a single global MAD over-cuts the shallow batch and under-cuts the deep one. Apply `is_outlier` within each `batch_key` group (the same per-batch logic the HVG step uses).

```python
flags = ['log1p_total_counts', 'log1p_n_genes_by_counts', 'pct_counts_in_top_20_genes']
adata.obs['outlier'] = adata.obs.groupby('sample', observed=True).apply(
    lambda g: (is_outlier(adata[g.index], flags[0], 5) | is_outlier(adata[g.index], flags[1], 5)
               | is_outlier(adata[g.index], flags[2], 5))).droplevel(0)
```

```r
# '^MT-' matches gene SYMBOLS; with Ensembl-ID feature names it matches nothing and the mito filter silently does nothing
seurat_obj[['percent.mt']] <- PercentageFeatureSet(seurat_obj, pattern = '^MT-')
VlnPlot(seurat_obj, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)
```

### QC Thresholds and Rationale

| Metric | Reference value | Rationale and caveat |
|--------|-----------------|----------------------|
| `min_genes` | 200 | Below this is mostly empty droplets / debris; raise for deep data |
| `log1p_total_counts` / `log1p_n_genes_by_counts` | 5 MAD | sc-best-practices loosens from scater's 3 MAD to avoid cutting real biology; filter on the log scale (depth is right-skewed) |
| `pct_counts_in_top_20_genes` | 5 MAD | High value flags low-complexity / dying cells |
| `pct_counts_mt` | 3 MAD plus hard >8% | Tissue-dependent: 5-20% typical, but cardiomyocytes/hepatocytes/muscle are constitutively high; nuclei are ~0-2% and any mito flags ambient |
| `min_cells` (genes) | 3 | Remove genes seen in too few cells to be informative |

Fixed cutoffs are a fast first pass for well-characterized tissue but silently delete valid populations; MAD-adaptive is the modern default; miQC (a mito-vs-detected-genes mixture model) helps when that relationship varies across samples.

### Mito and Dissociation Confounds

High mito is ambiguous: apoptosis co-occurs with low gene counts and apoptotic markers, while warm-dissociation stress co-occurs with immediate-early genes (FOS, JUN, JUNB, EGR1) and heat-shock proteins (HSPA1A/B) at normal gene counts. The IEG/HSP program creates a spurious "activated/stressed" cluster that passes every count/mito filter, is cell-type-specific in magnitude (so it does not cancel as a uniform batch effect), and overlaps real immune/stem activation, so naive removal can itself delete biology. Score the dissociation module per cell, then exclude those genes from HVG/clustering or flag and interpret cautiously; cold-protease digestion and single-nucleus assays reduce the artifact. For nuclei, standard mito thresholds are meaningless (baseline near zero) - lean on counts/genes outliers.

## Ambient RNA Removal

**Goal:** Remove cell-free "soup" mRNA that inflates off-target markers (hemoglobin everywhere in PBMCs, hepatocyte genes in non-hepatocytes) and fabricates co-expression.

**Approach:** Estimate the soup profile and a per-cell contamination fraction, then subtract; pick ONE tool and validate that a known-specific marker survives.

| Tool | Input | Needs empty droplets? | Strength | Fails / risk |
|------|-------|-----------------------|----------|--------------|
| SoupX (R) | Cell Ranger raw+filtered | Yes | Fast, interpretable rho, auto-estimate | `autoEstCont` fails on homogeneous data; single global soup wrong when ambient is heterogeneous |
| CellBender (Python, GPU) | RAW h5 | Yes (core of model) | Deep generative; removes ambient + barcode noise; also does cell-calling; strong on nuclei | Over-removes real low-abundance genes at high `--fpr`; black-box; slow |
| DecontX (R, celda) | Filtered cells | No | No raw needed; easy SCE/Seurat integration | Relies on cluster purity |

```r
library(SoupX)
sc <- load10X('cellranger_outs/')          # needs BOTH raw and filtered
sc <- autoEstCont(sc)                       # estimates contamination fraction rho
counts_adj <- adjustCounts(sc, roundToInt = TRUE)   # output is non-integer by default; round for NB models
```

SoupX and CellBender disagree on what "ambient" is: SoupX subtracts a per-cell scalar of a single global soup profile; CellBender learns a probabilistic per-droplet background in a generative model. There is no consensus on which is better - CellBender is more powerful and more dangerous. Subtracting a shared soup vector from every cell can manufacture artificial negative correlations and zero out genes cells genuinely lacked, so validate. Matters most for solid tumors, snRNA-seq, and blood. Do not stack tools; double-correction compounds over-removal.

## Normalization

**Goal:** Remove per-cell depth bias and stabilize variance so high-expression genes do not dominate PCA/kNN distances.

**Approach:** Default to shifted-log; reach for scran on shallow data and Pearson residuals on UMI count models; never normalize already-normalized data.

```python
adata.layers['counts'] = adata.X.copy()                  # stash raw before normalizing (HVG/doublets need it)
sc.pp.normalize_total(adata)                             # target_sum=None scales each cell to the dataset MEDIAN; pass target_sum=1e4 for the historical, arbitrary CP10k
sc.pp.log1p(adata)                                       # natural-log(1+x); the variance-stabilizing transform (no target_sum argument)
```

```r
seurat_obj <- NormalizeData(seurat_obj, normalization.method = 'LogNormalize', scale.factor = 10000)
# or variance-stabilized: seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
```

| Method | Model / assumption | Use when | Fails when |
|--------|--------------------|----------|------------|
| Shifted-log (CP10k / median) | Size-factor + log1p; constant total mRNA | General default; strong, fast, defensible | Composition-divergent types (plasma, cycling) distort fold-changes |
| scran deconvolution | Pooled size factors robust to composition | Low-depth, high-dropout, plate-based | R-only; needs pre-clustering (`quickCluster`); factors can go negative |
| sctransform v1/v2 | NB regularized regression (Pearson residuals) | Seurat depth removal for HVG/viz | Slow; off the count scale; v1 overfits theta (use v2) |
| Analytic Pearson residuals | `r=(x-mu)/sqrt(mu+mu^2/theta)` | UMI HVG+PCA without ad-hoc steps | Experimental; residual variance depends on theta and depth; clip to +/-sqrt(n) |

Ahlmann-Eltze and Huber 2023 found plain shifted-log + PCA performs as well as or better than sctransform, Pearson residuals, and GLM-PCA on kNN-overlap recovery, so shifted-log is the defensible default and the sophisticated methods are "use if preferred," not mandated. Because methods genuinely compete here, verify current best practice against the installed tool's docs before committing. Normalize raw counts exactly once and keep the transform consistent across HVG, scaling, and PCA.

## Highly Variable Genes

**Goal:** Restrict PCA/clustering to genes carrying biological signal.

**Approach:** Select a flavor, then feed it the input type it expects - the single most consequential gotcha is that dispersion flavors want log-normalized data while `seurat_v3` and Pearson want RAW COUNTS.

```python
# seurat_v3 reads raw counts from a layer and REQUIRES n_top_genes; needs the scikit-misc package
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3', layer='counts')
```

| Flavor | Function | Input | n_top_genes required? | Extra dependency |
|--------|----------|-------|-----------------------|------------------|
| `seurat` (default) | `sc.pp.highly_variable_genes` | log-normalized | No | - |
| `cell_ranger` | `sc.pp.highly_variable_genes` | log-normalized | No | - |
| `seurat_v3` | `sc.pp.highly_variable_genes` | RAW counts | Yes | scikit-misc |
| `pearson_residuals` | `sc.experimental.pp.highly_variable_genes` | RAW counts | recommended | - |

Running `seurat_v3` on logged values, or `seurat` on raw counts, runs silently and yields garbage HVGs. The field is shifting toward binomial-deviance and Pearson-residual feature selection on raw counts because dispersion HVGs are sensitive to the upstream normalization choice. Set `batch_key` to compute HVGs per batch and avoid batch-specific technical genes.

## Scaling and Regressing Out

**Goal:** Optionally equalize gene weight in PCA, and remove unwanted covariates - both now discouraged as reflexive defaults.

**Approach:** Prefer PCA on log-normalized HVG without scaling; regress out only a validated, non-confounded covariate.

```python
sc.pp.scale(adata, max_value=10)        # max_value default is None (no clipping); 10 is an explicit choice to cap z-scores
# sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])   # scanpy itself warns this overcorrects
```

"Always regress out mito and total_counts" is folklore: those covariates are confounded with real cell identity and state (cycling cells legitimately have more RNA), so regressing them erases biology and can collapse data into a blob. Modern normalization already stabilizes depth; address unwanted variation with integration (Harmony, scVI) rather than linear regression. Scaling inflates lowly-expressed noisy genes; sc-best-practices runs PCA on the normalized layer directly.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| HVGs look random; clustering is mush | `seurat_v3` fed log-normalized (or `seurat` fed raw) | Feed each flavor its required input; use `layer='counts'` for `seurat_v3` |
| An entire healthy cell type disappeared | Flat mito cutoff deleted high-mito parenchyma | Use MAD/tissue-aware thresholds; inspect what was removed |
| Values inflated ~2x after re-running normalization | Normalized already-normalized data | Normalize raw once; restore from `layers['counts']` |
| `ModuleNotFoundError: skmisc` | `seurat_v3` needs scikit-misc | `pip install scikit-misc` |
| QC metrics missing from `.obs` | `calculate_qc_metrics` `inplace` defaults to False | Pass `inplace=True` |
| Proliferation / activation signal vanished | Regressed out `total_counts` / cell-cycle confounded with biology | Do not reflexively regress; validate the covariate is not confounded |
| New "stressed/transitional" cluster | Warm-dissociation IEG/HSP artifact | Score the dissociation module; exclude those genes from HVG/clustering |
| Off-target markers everywhere (Hb, Ig) | Ambient RNA contamination | Run SoupX/CellBender/DecontX on the raw matrix before QC |
| Spike to ~2x counts deflated other genes | Compositional see-saw from a few dominant genes | Use `exclude_highly_expressed=True` or scran; report relative, not absolute, expression |
| Almost all cells filtered / tiny survivor count | MAD ~ 0 on a low-variance, tiny, or nuclei sample (>50% share a value), so `is_outlier` flags every non-median cell | Assert `n_obs > 0` and a sane survival fraction; fall back to fixed cutoffs when MAD is ~0 |
| Mito filter removes nothing (percent.mt all 0) | `'^MT-'` pattern matched against Ensembl-ID feature names | Use gene symbols, or match the mito Ensembl IDs / a mito gene list |
| Shallow batch over-filtered, deep batch under-filtered | Global MAD thresholds across samples of differing depth | Compute `is_outlier` per `batch_key`/sample group |
| Batch effects baked into QC/soup/doublet calls | Merged samples before QC, ambient, and doublet steps | Run steps 1-5 per sample, then merge |

## Related Skills

- single-cell/data-io - load the raw matrix before preprocessing
- single-cell/doublet-detection - per-sample doublet calling around the QC step
- single-cell/clustering - PCA, neighbors, and clustering after preprocessing
- single-cell/batch-integration - correct batch effects instead of regressing them out
- single-cell/markers-annotation - find markers after clustering
- differential-expression/deseq2-basics - pseudobulk DE across samples (avoids single-cell pseudo-replication)

## References

- Heumos L, Schaar AC, Lance C, et al. (2023) Best practices for single-cell analysis across modalities. Nature Reviews Genetics 24:550-572. DOI 10.1038/s41576-023-00586-w
- Ahlmann-Eltze C, Huber W (2023) Comparison of transformations for single-cell RNA-seq data. Nature Methods 20:665-672. DOI 10.1038/s41592-023-01814-1
- Lun ATL, Bach K, Marioni JC (2016) Pooling across cells to normalize single-cell RNA sequencing data (scran). Genome Biology 17:75. DOI 10.1186/s13059-016-0947-7
- Hafemeister C, Satija R (2019) Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression (sctransform). Genome Biology 20:296. DOI 10.1186/s13059-019-1874-1
- Choudhary S, Satija R (2022) Comparison and evaluation of statistical error models for scRNA-seq (sctransform v2). Genome Biology 23:27. DOI 10.1186/s13059-021-02584-9
- Lause J, Berens P, Kobak D (2021) Analytic Pearson residuals for normalization of single-cell RNA-seq UMI data. Genome Biology 22:258. DOI 10.1186/s13059-021-02451-7
- Vallejos CA, Risso D, Scialdone A, Dudoit S, Marioni JC (2017) Normalizing single-cell RNA sequencing data: challenges and opportunities. Nature Methods 14(6):565-571. DOI 10.1038/nmeth.4292
- Osorio D, Cai JJ (2021) Systematic determination of the mitochondrial proportion in human and mouse tissues for scRNA-seq quality control. Bioinformatics 37(7):963-967. DOI 10.1093/bioinformatics/btaa751
- Hippen AA, Falco MM, Weber LM, et al. (2021) miQC: An adaptive probabilistic framework for quality control of single-cell RNA-seq data. PLoS Computational Biology 17(8):e1009290. DOI 10.1371/journal.pcbi.1009290
- Young MD, Behjati S (2020) SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data. GigaScience 9(12):giaa151. DOI 10.1093/gigascience/giaa151
- Fleming SJ, Chaffin MD, Arduini A, et al. (2023) Unsupervised removal of systematic background noise (CellBender remove-background). Nature Methods 20:1323-1335. DOI 10.1038/s41592-023-01943-7
- van den Brink SC, Sage F, Vertesy A, et al. (2017) Single-cell sequencing reveals dissociation-induced gene expression in tissue subpopulations. Nature Methods 14(10):935-936. DOI 10.1038/nmeth.4437
- Squair JW, Gautier M, Kathe C, et al. (2021) Confronting false discoveries in single-cell differential expression. Nature Communications 12:5692. DOI 10.1038/s41467-021-25960-2
