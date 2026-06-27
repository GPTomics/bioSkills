---
name: bio-single-cell-cell-annotation
description: Automated reference-based cell type annotation for single-cell RNA-seq using CellTypist, SingleR, Azimuth, scANVI, and scmap to transfer labels from a reference. Use when annotating cell types from a reference atlas or pretrained model, transferring labels onto a query, assessing prediction confidence and rejection, or triaging whether an unexpected cluster is a novel type versus a doublet, low-quality, or batch artifact.
tool_type: mixed
primary_tool: CellTypist
---

## Version Compatibility

Reference examples tested with: scanpy 1.10+, Seurat 5.0+, celltypist 1.6+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Automated Reference-Based Cell Annotation

**"Annotate my cells from a reference"** -> Transfer labels from an annotated reference or pretrained model onto query cells, with a calibrated confidence/rejection step.
- Python: `celltypist.annotate()` (pretrained LR) or `scvi-tools` scANVI label transfer
- R: `SingleR()` (correlation to a reference) or `RunAzimuth()` (anchor-based mapping)

## Governing principle

An annotation is a HYPOTHESIS, not a measurement. Reference-based annotators are closed-world: every query cell is forced toward the nearest label the reference contains, so a genuinely novel state gets the nearest wrong label - often with high apparent confidence. Reproducible is not the same as correct; an automated label inherits the reference's annotation errors, granularity, and tissue/donor/disease scope, and propagates them at scale wearing the authority of "automated."
Markers are context-dependent. A marker gene is a conditional statement: "marker X = type Y" means "in this tissue, platform, and processing, X enriches in Y relative to these other cells." A marker in blood may be expressed broadly in tumor; "canonical confirmation" using the same markers that defined the type is circular (recovering the prior, not new evidence). Treat reference labels and marker catalogs (PanglaoDB, CellMarker) as priors to be triangulated, never ground truth.
Before naming a new cell type, triage the four-way confusion in decreasing frequency: (1) doublets - two cell types summed, co-expressing mutually exclusive lineage markers; (2) low-quality/dying - high mito %, low gene count, ambient-dominated; (3) batch/technical - the cluster maps to one sample/lane/chemistry; (4) ambient-RNA contamination (SoupX/CellBender). Only after excluding all four is "novel cell type" admissible. The field is littered with "novel populations" that were doublets or stress artifacts.

This skill covers automated reference transfer. Manual marker discovery and hand-labeling live in single-cell/markers-annotation; the two are complementary - automate a first pass, confirm with markers, reserve expert curation for the final label and ambiguous populations.

## Choosing an annotation method

| Method | Model | Reference | World | Use when | Fails when |
|--------|-------|-----------|-------|----------|------------|
| CellTypist | Logistic regression (pretrained) | Pretrained immune/cross-tissue models | Closed (+probability) | Immune/PBMC, fast first pass, no R needed | Input not log1p CP10K-normalized; query far from training distribution |
| SingleR | Spearman correlation to reference | celldex bulk or single-cell refs | Closed (+pruning) | Bulk reference available, R workflow, per-cell scoring | Strong platform/chemistry shift vs reference; forces nearest label |
| Azimuth | Supervised PCA + anchor mapping | Curated Seurat atlases (PBMC, lung...) | Closed (+mapping.score) | A curated Azimuth reference matches the tissue | No matching reference; locked to provided atlases |
| scANVI / scArches | Semi-supervised VAE | Annotated atlas + raw counts | Closed (+latent uncertainty) | Strong query batch vs reference; mapping onto a large atlas | Training cost/hyperparameters; raw counts required |
| scmap | Nearest reference centroid/cell | Single-cell reference | Open (explicit unassigned) | An explicit rejection category is needed | Coarser resolution; threshold tuning |
| LLM (GPTCelltype) | Prompted from top markers | None (uses marker list) | Open-ish | Fast hypothesis from a marker table | Hallucination, non-reproducible, never sees expression |

No method escapes the closed-world limit except by an explicit reject/unassigned bin. When methods compete, verify current best practice and reference availability against installed docs before committing.

## Normalization requirements (silent-failure risk)

| Tool | Required input | Wrong input symptom |
|------|----------------|---------------------|
| CellTypist | log1p-normalized to 10,000 counts/cell (CP10K) | Confident but degraded/wrong labels, no error |
| SingleR | log-normalized expression (`logcounts`) | Distorted correlations |
| scANVI/scArches | RAW counts in a layer | Model trains on the wrong likelihood |
| Azimuth | raw counts (SCTransform applied internally) | Mapping QC degrades |

## CellTypist (Python)

**Goal:** Transfer labels from a pretrained model with cluster-level smoothing and a probability for rejection.

**Approach:** Normalize the query to CP10K log1p (the model's expected input), run `annotate` with `majority_voting` to reassign each over-clustered subgroup to its dominant label, then keep a per-cell confidence for filtering.

```python
import scanpy as sc
import celltypist
from celltypist import models

adata = sc.read_h5ad('clustered.h5ad')
adata.X = adata.layers['counts'].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

models.download_models(model='Immune_All_Low.pkl')
predictions = celltypist.annotate(adata, model='Immune_All_Low.pkl', majority_voting=True)
adata = predictions.to_adata()

adata.obs['cell_type'] = adata.obs['majority_voting']
adata.obs['uncertain'] = adata.obs['conf_score'] < 0.5
```

## SingleR (R)

**Goal:** Assign each cell by correlation to a reference and prune low-confidence calls.

**Approach:** Score each cell's Spearman correlation to reference profiles (per-label score is the 0.8 quantile), assign the max, fine-tune, then prune cells whose delta (assigned-label score minus median) falls >3 MADs below the delta distribution.

```r
library(SingleR)
library(celldex)
library(SingleCellExperiment)

sce <- as.SingleCellExperiment(seurat_obj)
ref <- celldex::HumanPrimaryCellAtlasData()

pred <- SingleR(test = sce, ref = ref, labels = ref$label.main, de.method = 'classic', fine.tune = TRUE)
seurat_obj$SingleR <- pred$labels
seurat_obj$SingleR_pruned <- pred$pruned.labels

plotScoreHeatmap(pred)
plotDeltaDistribution(pred)
```

Use `de.method='classic'` for bulk references and `de.method='wilcox'` for single-cell references. Cells pruned to NA are the rejection set; inspect the delta distribution rather than trusting a hard score cutoff.

## Azimuth (R/Seurat)

**Goal:** Map a query onto a curated reference atlas and transfer hierarchical labels with a mapping score.

**Approach:** Project query cells onto the supervised reference embedding via anchors, transfer l1/l2/l3 labels, and gate by `mapping.score` and `prediction.score`.

```r
library(Seurat)
library(Azimuth)

seurat_obj <- RunAzimuth(seurat_obj, reference = 'pbmcref')
seurat_obj$azimuth <- seurat_obj$predicted.celltype.l2
seurat_obj$azimuth_low_conf <- seurat_obj$predicted.celltype.l2.score < 0.7
```

## Rejection thresholds (calibrate, do not port)

| Tool | Rejection signal | Default-ish |
|------|------------------|-------------|
| SingleR | delta + `pruneScores(nmads=3)` | 3 MADs below delta distribution |
| CellTypist | `conf_score` / `p_thres` | 0.5 |
| scmap | max similarity | < 0.7 unassigned |
| Azimuth/scANVI | mapping.score / latent uncertainty | inspect per dataset |

A hard universal probability cutoff is not principled across models - inspect the score distribution and calibrate per dataset.

## Triage an unexpected cluster (before claiming novelty)

**Goal:** Decide whether a poorly-mapped cluster is a novel type or an artifact.

**Approach:** A whole cluster scoring low (vs scattered low-confidence cells) suggests "not in reference"; rule out doublets, low-quality, batch, and ambient before annotating de novo.

```python
import numpy as np

cluster_conf = adata.obs.groupby('leiden')['conf_score'].median()
suspect = cluster_conf[cluster_conf < 0.5].index.tolist()

qc = adata.obs.groupby('leiden')[['pct_counts_mt', 'n_genes_by_counts', 'predicted_doublet']].mean()
print(qc.loc[suspect])
batch_purity = adata.obs.groupby('leiden')['sample'].agg(lambda s: s.value_counts(normalize=True).max())
print(batch_purity.loc[suspect])
```

High mito or low gene count flags low-quality; doublet rate or co-expressed exclusive lineages flags doublets; near-1 batch purity flags a technical artifact. Only a low-confidence, QC-clean, batch-mixed cluster with coherent de-novo markers is a novel-type candidate.

## Validate predictions with markers

**Goal:** Confirm transferred labels against canonical markers (triangulation, not proof).

**Approach:** Dot-plot lineage markers grouped by predicted label and check the expected on/off pattern; disagreement between automated calls and markers flags cells to re-examine.

```r
canonical <- c('CD3D', 'CD8A', 'MS4A1', 'CD14', 'FCGR3A', 'NKG7', 'FCER1A')
DotPlot(seurat_obj, features = canonical, group.by = 'SingleR') + Seurat::RotatedAxis()
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Confident labels that contradict canonical markers | Closed-world: novel/absent state forced to nearest label | Add a reject bin; annotate de novo; do not trust labels outside the reference's domain |
| CellTypist labels degrade silently | Query not CP10K log1p normalized | Normalize to `target_sum=1e4` then `log1p` before annotate |
| CellTypist returns confident but nonsensical labels | Gene-ID space mismatch (query `var_names` are Ensembl IDs vs symbol-based model); few genes matched | Set `var_names` to gene symbols; check the matched-gene fraction reported by annotate before trusting labels |
| Reference labels look wrong everywhere | Platform/chemistry shift vs reference (domain shift) | Use a batch-modeling mapper (scANVI/scArches) or a matched reference |
| "Novel cell type" turns out artifactual | Doublet / low-quality / batch / ambient not excluded | Run the four-way triage before claiming novelty |
| Fine labels (CD4 Tcm vs Tem) unstable | Granularity finer than data or reference supports | Annotate hierarchically; report coarse labels confidently, fine as hypotheses |
| Two tools disagree on the same cells | Different references/granularity | Report consensus + flag disagreements as ambiguous; curate manually |

## Related Skills

- markers-annotation - Manual marker discovery and hand-labeling that complements automated transfer
- clustering - Cluster cells before annotating
- preprocessing - Normalize correctly for each annotator's expected input
- batch-integration - Reference mapping vs de-novo integration; closed-world caveats
- differential-abundance - Test whether annotated cell-type proportions changed between conditions
- pathway-analysis/go-enrichment - Functionally characterize a de-novo / novel population

## References

- Aran et al. 2019, Nat Immunol 20:163-172 - SingleR correlation-based reference annotation with delta-based pruning.
- Dominguez Conde et al. 2022, Science 376:eabl5197 - CellTypist logistic-regression cross-tissue immune annotation.
- Hao et al. 2021, Cell 184(13):3573-3587 - Azimuth / weighted-NN reference mapping and label transfer.
- Xu et al. 2021, Mol Syst Biol 17(1):e9620 - scANVI semi-supervised annotation with calibrated uncertainty.
- Kiselev, Yiu & Hemberg 2018, Nat Methods 15:359-362 - scmap projection with an explicit unassigned category.
- Hou & Ji 2024, Nat Methods 21(8):1462-1465 - GPT-4 / GPTCelltype marker-based annotation and its hallucination/reproducibility caveats.
