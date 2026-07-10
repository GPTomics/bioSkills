---
name: bio-workflows-cytometry-pipeline
description: End-to-end flow, spectral, and mass cytometry (CyTOF) pipeline from raw FCS files to differentially abundant/expressed cell populations. Orchestrates the read -> compensate/unmix -> transform -> QC -> doublet-removal -> cluster-or-gate -> annotate -> diffcyt DA/DS chain with flowCore/CATALYST/diffcyt, branching on instrument type and on clustering-vs-gating. Use when processing a cytometry experiment end-to-end, deciding the pipeline path for an instrument, or wiring the flow-cytometry component skills into one analysis with valid sample-level statistics.
tool_type: r
primary_tool: CATALYST
workflow: true
depends_on:
  - flow-cytometry/fcs-handling
  - flow-cytometry/compensation-transformation
  - flow-cytometry/cytometry-qc
  - flow-cytometry/doublet-detection
  - flow-cytometry/bead-normalization
  - flow-cytometry/gating-analysis
  - flow-cytometry/clustering-phenotyping
  - flow-cytometry/differential-analysis
qc_checkpoints:
  - after_load: "Read RAW (transformation=FALSE); >~10K cells/sample for stable per-sample frequencies"
  - after_compensation: "Compensate/unmix on LINEAR data before transform; single-stain controls >= as bright as sample; FMO for boundaries"
  - after_qc: "Margins removed BEFORE density QC; doublets removed BEFORE clustering; dead cells gated"
  - after_clustering: "10-30 metaclusters (over-provision then merge); cluster on TYPE markers only, test STATE in DS"
  - after_testing: ">=2-3 biological replicates/group (the sample is the unit); batch modeled in the design; BH FDR across clusters (and clusters x markers for DS)"
---

## Version Compatibility

Reference examples tested with: CATALYST 1.26+, diffcyt 1.22+, FlowSOM 2.10+, flowCore 2.14+, flowWorkspace 4.14+, flowStats 4.14+, edgeR 4.0+, limma 3.58+, ggplot2 3.5+; Python (partial alt) flowkit 1.1+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt rather than retrying. Each stage defers depth to its component skill.

# Flow Cytometry Pipeline

**"Process my cytometry data from FCS to differential populations"** -> read raw -> compensate/unmix -> transform -> QC -> remove doublets -> cluster (or gate) -> annotate -> test DA/DS, with the sample as the unit of inference.
- R: `flowCore` + `CATALYST::prepData/cluster/runDR` + `diffcyt::diffcyt()`

## The Single Most Important Modern Insight -- A Pipeline Is a Chain of Irreversible Decisions, and the Unit of Inference Is the Sample

Each early choice silently gates the validity of the final test: reading raw (not log-linearized), compensating BEFORE transforming, removing margin events before density QC, assigning type-vs-state markers correctly, and removing doublets before clustering. None of these is recoverable downstream - a doublet clustered as a "double-positive," a state marker used for clustering, or an uncompensated channel becomes a false population that the differential test then "confirms." The second critical thread is that the SAMPLE/subject, not the cell, is the experimental unit: diffcyt aggregates cells to per-sample-per-cluster counts (DA) and medians (DS) before testing, so biological replication (>= 2-3 per group) is mandatory and a per-cell test is invalid. Two normalization layers sit at different points in the pipeline - EQ-bead drift correction on raw counts at the very front (CyTOF), and CytoNorm cross-batch harmonization on transformed data before the analytical clustering (its internal FlowSOM clustering is part of the batch model, not the analysis) - and conflating them is a classic error.

## Decision Tree: Which Path

| Situation | Path | Why |
|-----------|------|-----|
| Conventional fluorescence flow | compensate (`$SPILLOVER`/flowStats) -> logicle -> ... | optical spillover; logicle handles negatives |
| Spectral cytometer (Aurora/ID7000) | UNMIX (not compensate) -> arcsinh ~150 | overdetermined system; fluorescence-scale |
| Mass cytometry (CyTOF) | EQ-bead normalize (raw) -> arcsinh cofactor 5 -> `compCytof` if needed | metals barely spill (~1-4%); drift correction first |
| High-dim discovery, no prior gates | cluster (FlowSOM via CATALYST) | scales; finds unexpected populations |
| Well-defined populations / rare events (MRD) | hierarchical gating (openCyto) | interpretable; clustering fails for ultra-rare |
| Multi-batch / multi-day | anchor sample per batch -> CytoNorm (normalize transformed data before analytical clustering) | model batch in the design for inference |

## Pipeline Overview

```
FCS -> compensate/unmix -> transform -> QC (margins, time, dead) -> doublets
     -> [ cluster (FlowSOM) | gate (openCyto) ] -> annotate -> diffcyt DA/DS -> report
EQ-bead drift normalization (CyTOF) runs on raw counts BEFORE everything; CytoNorm runs on transformed data and its normalized output feeds the cluster/gate step.
```

## 1. Panel, Metadata, and Load

**Goal:** Define the type/state panel and sample metadata, then load FCS.

**Approach:** Panel `marker_class` drives everything downstream (type clusters, state is tested); metadata keys samples to condition/subject. See flow-cytometry/fcs-handling.

```r
library(CATALYST); library(diffcyt); library(flowCore); library(ggplot2)

panel <- data.frame(
  fcs_colname = c('FSC-A','SSC-A','CD45','CD3','CD4','CD8','CD19','CD14','Ki67','IFNg'),
  antigen     = c('FSC','SSC','CD45','CD3','CD4','CD8','CD19','CD14','Ki67','IFNg'),
  marker_class = c('none','none','type','type','type','type','type','type','state','state'))
md <- data.frame(file_name = list.files('data', pattern = '\\.fcs$'),
                 sample_id = paste0('S', 1:8),
                 condition = rep(c('Control','Treatment'), each = 4),
                 patient_id = rep(paste0('P', 1:4), 2))
fs <- read.flowSet(file.path('data', md$file_name), transformation = FALSE, truncate_max_range = FALSE)
```

## 2. Compensate / Unmix, then Transform

**Goal:** Remove spillover on linear data, then variance-stabilize.

**Approach:** Conventional flow compensates (matrix before transform); CyTOF skips fluorescence compensation and uses cofactor 5; spectral unmixes then uses ~150. See flow-cytometry/compensation-transformation.

```r
spill <- spillover(fs[[1]]); spill <- spill[[which(!vapply(spill, is.null, logical(1)))[1]]]  # first POPULATED matrix; FACS stores it under SPILL/$SPILLOVER, not always [[1]]
fs_comp <- compensate(fs, spill)                        # conventional flow; CyTOF: omit or use compCytof
COFACTOR <- 150                                          # 5 for CyTOF, ~150 for fluorescence/spectral
sce <- prepData(fs_comp, panel, md, transform = TRUE, cofactor = COFACTOR, FACS = TRUE)
```

## 3. QC (order matters)

**Goal:** Remove margin/boundary events and time anomalies before any density step.

**Approach:** Margins first, then time-based cleaning; on CyTOF, EQ-bead drift correction happens upstream on raw counts. See flow-cytometry/cytometry-qc and flow-cytometry/bead-normalization.

```r
# per-sample sanity + sample-similarity MDS (flag outlier samples)
plotExprs(sce, color_by = 'condition'); pbMDS(sce, color_by = 'condition')
# event-level cleaning runs per-FCS upstream: PeacoQC::RemoveMargins() -> PeacoQC()/flowAI on transformed data
```

## 4. Remove Doublets

**Goal:** Drop aggregates before clustering so they don't form phantom double-positives.

**Approach:** Flow uses the FSC-A vs FSC-H diagonal; CyTOF uses DNA intercalator + Gaussian/Event_length. See flow-cytometry/doublet-detection.

```r
# CyTOF (FACS=TRUE retained Event_length on the arcsinh scale):
e <- assay(sce, 'exprs')
if (all(c('DNA1','Event_length') %in% rownames(sce))) {
  keep <- e['DNA1', ] > quantile(e['DNA1', ], 0.05) &
          e['Event_length', ] <= quantile(e['Event_length', ], 0.99)
  sce <- sce[, keep]
}
```

## 5. Cluster (FlowSOM) or Gate

**Goal:** Define populations by unsupervised clustering on TYPE markers (discovery) or hierarchical gating (defined/rare).

**Approach:** `cluster()` wraps FlowSOM+ConsensusClusterPlus; over-provision the grid, set a seed. See flow-cytometry/clustering-phenotyping (clustering) and flow-cytometry/gating-analysis (gating).

```r
sce <- cluster(sce, features = 'type', xdim = 10, ydim = 10, maxK = 20, seed = 42)
```

## 6. Annotate and Visualize Structure

**Goal:** Label metaclusters from marker medians; embed for display only.

**Approach:** Median heatmap drives annotation; UMAP colors by cluster but is never used to define or quantify populations.

```r
plotExprHeatmap(sce, features = 'type', by = 'cluster_id', k = 'meta20', scale = 'last')
sce <- runDR(sce, dr = 'UMAP', features = 'type', cells = 2000)
plotDR(sce, dr = 'UMAP', color_by = 'meta20')
```

## 7. Differential Abundance and State

**Goal:** Test which populations change in frequency (DA) or state-marker expression (DS) between conditions.

**Approach:** The `diffcyt()` wrapper aggregates to the sample level; results live in `res$res`. See flow-cytometry/differential-analysis.

```r
design   <- createDesignMatrix(ei(sce), cols_design = 'condition')
contrast <- createContrast(c(0, 1))                        # Treatment vs Control
res_DA <- diffcyt(sce, clustering_to_use = 'meta20', analysis_type = 'DA',
                  method_DA = 'diffcyt-DA-edgeR', design = design, contrast = contrast)
res_DS <- diffcyt(sce, clustering_to_use = 'meta20', analysis_type = 'DS',
                  method_DS = 'diffcyt-DS-limma', design = design, contrast = contrast)
da <- as.data.frame(SummarizedExperiment::rowData(res_DA$res))   # cluster_id, logFC, p_val, p_adj
```

## 8. Visualize Results and Export

**Goal:** Summarize significant populations and persist results.

**Approach:** Pass the inner result object (`res$res`) to plotting; export tables and the SCE.

```r
plotDiffHeatmap(sce, res_DA$res, all = TRUE, fdr = 0.05)
plotAbundances(sce, k = 'meta20', by = 'cluster_id', group_by = 'condition')
write.csv(da, 'da_results.csv', row.names = FALSE); saveRDS(sce, 'cytometry_analysis.rds')
```

## Paired / Repeated-Measures Variant

**Goal:** Account for within-subject correlation (pre/post on the same donor).

**Approach:** Use a GLMM with a random effect for subject (NOT voom, which is fixed-effects only).

```r
formula <- createFormula(ei(sce), cols_fixed = 'condition', cols_random = 'patient_id')
res_DA <- diffcyt(sce, clustering_to_use = 'meta20', analysis_type = 'DA',
                  method_DA = 'diffcyt-DA-GLMM', formula = formula, contrast = createContrast(c(0, 1)))
```

## Manual Gating Path (alternative to clustering)

**Goal:** Define populations by a reproducible hierarchy when they are well-defined or rare.

**Approach:** Build a GatingSet on transformed data; recompute after adding gates. See flow-cytometry/gating-analysis.

```r
library(flowWorkspace)
tl <- estimateLogicle(fs_comp[[1]], colnames(spill))
gs <- GatingSet(transform(fs_comp, tl))
# add openCyto template or manual gates (time -> debris -> singlets -> live -> lineage), then:
recompute(gs); gs_pop_get_stats(gs, type = 'count')
```

## Python Alternative (FlowKit) -- partial

**Goal:** Read, compensate, and gate in Python where an R pipeline is not an option.

**Approach:** FlowKit covers IO/compensation/GatingML; there is NO Python equivalent for diffcyt DA/DS, so the differential step stays in R (or bridge via readfcs -> AnnData -> scanpy for clustering only).

```python
import flowkit as fk
sample = fk.Sample('sample.fcs')
sample.apply_compensation(sample.metadata['spillover'])    # FlowKit lowercases + strips $ from keys, so $SPILLOVER -> 'spillover' (not 'spill'); use FlowKit's API, not a hand-rolled inverse
df = sample.as_dataframe(source='comp')
```

## Per-Stage Failure Modes

### Per-cell pseudoreplication
**Trigger:** testing across all cells. **Mechanism:** cells are not independent replicates. **Symptom:** p ~ 1e-40 from few subjects. **Fix:** diffcyt aggregates to sample level; require >= 2-3 replicates/group.

### Clustering on state markers
**Trigger:** activation/phospho markers in `features`. **Mechanism:** state contaminates lineage identity. **Symptom:** activated/resting splits of one type. **Fix:** cluster on `type`; test state in DS.

### Doublets / wrong cofactor / uncompensated input
**Trigger:** skipping doublet removal, cofactor 5 on fluorescence, or clustering raw data. **Mechanism:** phantom double-positives, compressed dim markers, spillover-dominated distances. **Symptom:** non-reproducible "novel" populations. **Fix:** remove doublets first; cofactor 5 (CyTOF) / 150 (fluorescence); compensate+transform before clustering.

### Batch cleaned instead of modeled
**Trigger:** CytoNorm-ing then testing naively, or batch confounded with condition. **Mechanism:** over-correction / non-identifiability. **Symptom:** attenuated or fabricated effects. **Fix:** model batch in the design; if batch == condition, no rescue.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| arcsinh cofactor 5 (CyTOF) / ~150 (fluorescence) | Nowicka 2017 *F1000Res* 6:748 | matches platform noise scale |
| >= 2-3 biological replicates per group | Weber 2019 *Commun Biol* 2:183 | minimum for a valid DA/DS error term |
| > ~10K cells per sample | community | stable per-sample cluster frequencies |
| 10-30 metaclusters typical (maxK=20 default) | Weber & Robinson 2016 *Cytometry A* 89:1084 | over-provision then merge |
| BH FDR across clusters (and clusters x markers for DS) | diffcyt | many simultaneous tests |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `testDA_edgeR(sce, ...)` not found / wrong | fabricated signature | use the `diffcyt()` wrapper; results in `res$res` |
| `compensate()` errors / silent NULL | `spillover(ff)` returns a 3-slot list; the matrix is often under `SPILL`/`$SPILLOVER`, not `[[1]]` | select the first non-null slot, not positional `[[1]]` |
| empty DS results | state markers not flagged | set `marker_class='state'` in the panel |
| paired design ignored | used fixed-effect method | `diffcyt-DA-GLMM` with a random effect |

## References

- Weber 2019 *Commun Biol* 2:183 — diffcyt DA/DS framework.
- Nowicka 2017 *F1000Research* 6:748 — CATALYST CyTOF workflow; type/state, cofactor 5.
- Weber & Robinson 2016 *Cytometry A* 89(12):1084-1096 — FlowSOM clustering benchmark.
- Van Gassen 2020 *Cytometry A* 97(3):268-278 — CytoNorm cross-batch normalization.
- Hurlbert 1984 *Ecol Monogr* 54(2):187-211 — pseudoreplication (sample is the unit).

## Related Skills

- flow-cytometry/fcs-handling - Read FCS and map channels
- flow-cytometry/compensation-transformation - Compensate/unmix and transform
- flow-cytometry/cytometry-qc - Time/margin/dead-cell QC
- flow-cytometry/doublet-detection - Singlet discrimination
- flow-cytometry/bead-normalization - EQ-bead drift and CytoNorm batch correction
- flow-cytometry/gating-analysis - Hierarchical/automated gating path
- flow-cytometry/clustering-phenotyping - FlowSOM clustering and annotation
- flow-cytometry/differential-analysis - diffcyt DA/DS testing
- single-cell/clustering - Related graph-clustering for scRNA-seq
