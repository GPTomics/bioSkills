---
name: bio-single-cell-markers-annotation
description: Detect cluster marker genes and assign manual cell type labels in single-cell RNA-seq using Scanpy (Python) and Seurat (R). Use when finding genes that distinguish clusters, ranking markers for annotation, scoring gene signatures, hand-labeling clusters, or deciding between Wilcoxon marker ranking and pseudobulk condition DE.
tool_type: mixed
primary_tool: Seurat
---

## Version Compatibility

Reference examples tested with: scanpy 1.10+, Seurat 5.0+, anndata 0.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Marker Gene Detection and Manual Annotation

**"Find marker genes for my clusters"** -> Rank genes that separate each cluster from the rest, then map clusters to cell types using canonical markers.
- Python: `sc.tl.rank_genes_groups()` -> filter by effect size + fraction expressing -> `adata.obs[...].map(labels)`
- R: `Seurat::FindAllMarkers()` -> filter by `avg_log2FC` + `pct.1`/`pct.2` -> `RenameIdents()`

## Governing principle

Marker detection is descriptive ranking, NOT inference. Two distinct questions get sloppily called "DE" and must never be conflated: (1) marker detection - "which genes are higher in cluster X vs the rest?" is a ranking/annotation task where the cell is the unit and Wilcoxon is an acceptable heuristic; (2) condition DE - "which genes change in cell type X between treated and control?" is a population claim that requires biological replicates, where the unit of replication is the sample/donor, not the cell. Question 2 must use pseudobulk (aggregate raw counts per sample x cell type, then DESeq2/edgeR/limma-voom); treating cells as replicates is pseudoreplication and inflates false positives by orders of magnitude (Squair 2021).
Post-clustering marker p-values are double-dipping. Clusters were defined to maximize between-group separation, so testing those same clusters for markers tests a hypothesis built from the data used to test it. The Wilcoxon/t-test null assumes fixed a-priori labels; under a single homogeneous population the statistic does not follow its nominal null, type-I error approaches 1 as resolution rises, and BH correction does nothing because the p-values are invalid before correction. Cluster-marker p-values are descriptive labels, never evidence that a cluster is a real cell type. Rank and filter markers by effect size and fraction-expressing, not by p-value; a gene can be "significant" at p=1e-40 (n is thousands) yet useless as a marker (60% in-group vs 55% out-group).

This skill covers marker discovery for clusters plus manual labeling. Automated reference-based label transfer (SingleR, CellTypist, Azimuth, scANVI) lives in single-cell/cell-annotation. Cross-condition compositional change lives in single-cell/differential-abundance.

## Choosing a marker / DE method

| Method | Question answered | Use when | Fails when |
|--------|-------------------|----------|------------|
| Wilcoxon rank-sum (presto) | Rank cluster markers | Default for labeling a cluster vs rest; fast, non-parametric | Quoted as inference; double-dipping on the clustered data |
| t-test | Rank cluster markers | Quick first pass; scanpy `method=None` default | Heavy-tailed sparse counts violate normality; less robust than Wilcoxon |
| logistic regression (`logreg`/`LR`) | Markers controlling covariates | Need to adjust for batch/covariate when ranking | Slow; needs enough cells; still descriptive |
| ROC (`roc`, Seurat) | Classification power per gene | Want an AUC ranking of marker discriminativeness | No p-value; pure ranking |
| ClusterDE / count splitting | Are the cluster's markers real (FDR-honest)? | Validating that a split is not spurious before naming it | Adds a synthetic-null / data-thinning step; assumptions on the noise model |
| Pseudobulk + DESeq2/edgeR/limma-voom | Condition DE within a cell type | Treatment vs control with >=3 biological replicates per condition | n=1/condition (dispersion unidentifiable); cells-as-replicates |

Marker tools (`rank_genes_groups`, `FindMarkers`) will technically run a treatment-vs-control contrast cell-by-cell and return tidy tiny p-values. That is statistically invalid for a population claim. The tool not stopping the user is why this error is so common. When methods compete, verify current defaults against installed docs.

## Defaults that bite (verify before trusting tutorials)

| Tool | Folklore | Actual default |
|------|----------|----------------|
| scanpy `rank_genes_groups` | Defaults to Wilcoxon | `method=None` resolves to `t-test`; pass `method='wilcoxon'` explicitly |
| Seurat v5 `FindMarkers` `logfc.threshold` | 0.25 | 0.1 in v5 (was 0.25 in v4); permissive, returns more hits |
| Seurat v5 `FindMarkers` `min.pct` | 0.1 | 0.01 in v5 (was 0.1 in v4) |
| Seurat `test.use='wilcox'` | Always fast | Fast only if `presto` is installed; else silent slow base-R fallback |
| Pseudobulk input | Normalized/log values | Aggregate RAW counts (summed), never normalized |

## Scanpy marker detection

**Goal:** Rank cluster-specific markers and filter them by specificity, not p-value alone.

**Approach:** Run Wilcoxon explicitly (scanpy's default is t-test), pull results to a DataFrame with `pts=True` for in/out fraction, then keep genes with a large positive log fold change and a high in-group / low out-group fraction.

```python
import scanpy as sc

adata = sc.read_h5ad('clustered.h5ad')

sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', pts=True, corr_method='benjamini-hochberg')
markers = sc.get.rank_genes_groups_df(adata, group=None)

specific = markers[(markers['logfoldchanges'] > 1) & (markers['pct_nz_group'] > 0.5) & (markers['pct_nz_reference'] < 0.25)]
print(specific.groupby('group').head(10)[['group', 'names', 'logfoldchanges', 'pct_nz_group', 'pct_nz_reference']])
```

## Seurat marker detection

**Goal:** Rank markers per cluster and keep specific ones for labeling.

**Approach:** Run `FindAllMarkers` with `only.pos=TRUE`, install presto so Wilcoxon is fast, then rank within cluster by `avg_log2FC` and require a `pct.1`-`pct.2` gap.

```r
library(Seurat)
library(dplyr)

all_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1)

specific <- all_markers %>%
    filter(p_val_adj < 0.05, avg_log2FC > 1, (pct.1 - pct.2) > 0.2) %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)
print(specific)
```

Seurat v5 lowers thresholds to 0.1/0.01, so explicit `logfc.threshold=0.25` and a `pct.1-pct.2` filter restore a marker-grade (specific) shortlist from a permissive run.

## Gene signature scoring

**Goal:** Score each cell for a curated panel without library-size confounding.

**Approach:** Both tools subtract an expression-binned control set; thresholds are dataset-relative and must never be ported as absolute cutoffs.

```python
t_cell_panel = ['CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B']
sc.tl.score_genes(adata, gene_list=t_cell_panel, ctrl_size=50, n_bins=25, score_name='T_cell_score')
```

```r
seurat_obj <- AddModuleScore(seurat_obj, features = list(c('CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B')), ctrl = 100, name = 'T_cell_score')
```

scanpy uses 25 control bins, Seurat uses 24 by default (both follow Tirosh 2016) - a real cross-ecosystem non-reproducibility source for small panels.

## Cell-cycle scoring

**Goal:** Assign each cell an S and G2/M score and a phase, to diagnose (and optionally regress) cell-cycle-driven structure.

**Approach:** Score the Tirosh S and G2/M gene panels; both tools ship the lists. Regression is optional and confounded with biology (cycling is a real state in proliferating populations) - diagnose first and regress only when the cycle is a confound, not reflexively.

```python
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
```

```r
seurat_obj <- CellCycleScoring(seurat_obj, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
```

Provide `s_genes`/`g2m_genes` as the Tirosh 2016 panels (Seurat's `cc.genes.updated.2019` exposes both lists directly); scanpy ships no built-in list, so load the panels from the reference or a regev-lab gene file.

## Manual cluster labeling

**Goal:** Map cluster ids to cell type names after marker inspection.

**Approach:** Build a cluster->label dictionary from canonical-marker evidence, map it onto cells, and flag unmapped clusters rather than silently dropping them.

```python
cluster_labels = {'0': 'CD4 T', '1': 'CD14 Mono', '2': 'B', '3': 'CD8 T', '4': 'NK', '5': 'FCGR3A Mono'}
adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_labels).fillna('Unassigned')
```

```r
new_ids <- c('0' = 'CD4 T', '1' = 'CD14 Mono', '2' = 'B', '3' = 'CD8 T', '4' = 'NK', '5' = 'FCGR3A Mono')
seurat_obj <- RenameIdents(seurat_obj, new_ids)
seurat_obj$cell_type <- Idents(seurat_obj)
```

## Condition DE the correct way (pseudobulk)

**Goal:** Test which genes change between conditions within a cell type, with valid FDR.

**Approach:** Aggregate RAW counts to one profile per sample x cell type, then hand the count matrix to a bulk engine (DESeq2/edgeR/limma-voom) which estimates dispersion across biological replicates. Run each cell type separately so a one-cell-type effect is not diluted.

```python
import scanpy as sc

cell_type = adata[adata.obs['cell_type'] == 'CD14 Mono']
pseudobulk = sc.get.aggregate(cell_type, by='sample', func='sum')
counts_df = pseudobulk.layers['sum']
```

```r
pb <- AggregateExpression(seurat_obj, group.by = c('cell_type', 'sample'), assays = 'RNA', layer = 'counts')$RNA
```

Pull the summed counts slot, build a sample-level design (condition + covariates), and run DESeq2/edgeR; see differential-expression/deseq2-basics for the modeling step. Never run DE on batch-corrected or normalized expression.

## Canonical PBMC markers (context-dependent, validate per dataset)

| Cell type | Markers | Cell type | Markers |
|-----------|---------|-----------|---------|
| CD4 T | CD3D, CD4, IL7R | NK | NKG7, GNLY, NCAM1 |
| CD8 T | CD3D, CD8A, CD8B | CD14 Mono | CD14, LYZ, S100A8 |
| B | MS4A1, CD79A, CD19 | FCGR3A Mono | FCGR3A, MS4A7 |
| DC | FCER1A, CST3 | Platelet | PPBP, PF4 |

A marker is a conditional statement, not a property of a gene: a marker in blood may be expressed broadly in tumor, and "vs rest" markers depend on what "rest" is. Re-validate any ported panel.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Thousands of "significant" markers between two visually-similar clusters | Over-clustering + double-dipping inflation | Significance-test the split (scSHC/ClusterDE) or merge; never quote raw marker p-values as proof of a cell type |
| Marker p-values used as evidence clusters are real | Selective-inference violation; BH cannot fix invalid p-values | Report markers as descriptive labels; validate identity with orthogonal markers |
| Condition DE returns huge gene lists, none replicate | Cells treated as replicates (pseudoreplication) | Aggregate to pseudobulk per sample x cell type; test across donors |
| `FindAllMarkers` hangs for minutes | presto not installed; slow base-R Wilcoxon | `install.packages('presto')` (or `remotes::install_github('immunogenomics/presto')`) |
| Same top markers in every cluster | Resolution too high; clusters split one population | Lower resolution / merge; check stability |
| Gene cutoff ported from another dataset misclassifies cells | Module scores are dataset-relative | Set thresholds from this dataset's score distribution |
| NaN / degenerate logFC and p-values from marker ranking | Only one cluster present, so the "vs rest" reference is empty | Marker ranking needs >=2 groups; subcluster the population or report it as a single homogeneous type |
| "DE genes" between conditions but no gene changed per cell | Subpopulation proportions shifted (compositional confound) | Pair condition DE with single-cell/differential-abundance |

## Related Skills

- clustering - Cluster cells before finding markers
- preprocessing - Normalize and select features before marker detection
- cell-annotation - Automated reference-based label transfer (complements manual marker labeling)
- differential-abundance - Test whether cell-type proportions changed between conditions
- differential-expression/deseq2-basics - Pseudobulk condition DE engine for the aggregated counts
- differential-expression/de-results - Shrink, filter, and interpret pseudobulk DE results
- pathway-analysis/go-enrichment - Functional interpretation of marker / DE gene lists

## References

- Squair et al. 2021, Nat Commun 12:5692 - cells-as-replicates inflate false positives; top DE methods aggregate to pseudobulk.
- Crowell et al. 2020, Nat Commun 11:6077 - muscat; pseudobulk gives well-calibrated FDR for multi-sample multi-condition DS analysis.
- Neufeld et al. 2024, Biostatistics 25(1):270-287 - count splitting / valid post-clustering inference; Poisson thinning breaks under overdispersion.
- Lee & Han 2024, Bioinformatics 40(8):btae498 - properly-offset pseudobulk is statistically equivalent to a GLMM.
- Tirosh et al. 2016, Science 352:189-196 - control-set module scoring underlying score_genes / AddModuleScore.
