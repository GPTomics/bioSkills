---
name: bio-machine-learning-atlas-mapping
description: Maps query single-cell data onto reference atlases and transfers cell-type labels using scArches surgery (scVI/scANVI), Symphony, Azimuth, CellTypist, scPoli, popV, and foundation models, with explicit out-of-distribution and label-transfer uncertainty. Use when annotating new single-cell datasets against a pre-trained reference, deciding which mapping method fits, or judging whether transferred labels are trustworthy. For de novo clustering and manual annotation see single-cell/cell-annotation; for batch integration without a reference see single-cell/batch-integration.
tool_type: python
primary_tool: scvi-tools
---

## Version Compatibility

Reference examples tested with: anndata 0.10+, scanpy 1.10+, scvi-tools 1.1+, scikit-learn 1.3+, celltypist 1.6+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

scvi-tools 1.x has had API churn around minified/registered models and the semi-supervised loader. Confirm `scvi.__version__` and `help(scvi.model.SCANVI.load_query_data)` before relying on argument names. If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Reference Mapping and Label Transfer for Single-Cell Data

**"Annotate my scRNA-seq query against a reference atlas"** -> Project query cells into a fixed reference latent space and transfer labels, then gate every label on an out-of-distribution signal.
- scArches surgery: `scvi.model.SCANVI.prepare_query_anndata()` -> `load_query_data()` -> `train(plan_kwargs={'weight_decay': 0.0})` -> `predict(soft=True)`
- Linear/fixed reference: Symphony (R) or Azimuth/Seurat anchor transfer
- Pure classifier (no embedding): CellTypist logistic regression

## The Single Most Important Modern Insight -- Mapping Is a Projection, Not an Annotation Oracle

Reference mapping indicates *where in a fixed reference manifold a query cell lands*; it does not establish whether that location is biologically meaningful for the query. The projection always succeeds geometrically: every query cell is assigned its nearest reference label whether or not that label is true. A softmax over reference classes is normalized to sum to 1 and therefore *cannot* express "I have never seen this cell" -- a hepatocyte handed to an immune reference gets confidently called a T cell. Every failure mode below is a corollary of confusing projection with annotation.

The operational consequence: a transferred label is untrustworthy until paired with an out-of-distribution / label-transfer-uncertainty signal that measures *whether the cell belongs to the reference at all*, which is a different quantity from the prediction probability that measures *which reference label*. Conflating these two is the field's most common error.

## Methods Taxonomy

| Method | Model class | Needs reference model? | Novel-cell-type signal | Best when | Fails when |
|--------|-------------|------------------------|------------------------|-----------|------------|
| scVI + scArches surgery | Conditional VAE; surgery freezes reference weights, fits query-batch nodes | Yes (saved scVI model) | None intrinsic; add kNN uncertainty / OOD distance | Unseen batch; want a de novo joint embedding then cluster query yourself | Expected to *label* (it only embeds); reference lacks query biology |
| scANVI + scArches surgery | Semi-supervised VAE (scVI latent + label classifier head) | Yes (scANVI, often `from_scvi_model`) | Classifier softmax (overconfident); kNN-on-latent uncertainty | Reference well-labeled, query is the *same* tissue/biology | Semi-supervised leakage carves latent; novel states confidently mislabeled |
| Symphony | Linear Harmony soft-cluster mixture; query projected into fixed reference | Yes (compressed reference object) | Per-cell Mahalanobis distance to soft-cluster centroids | Seconds-scale, deterministic, CPU-only, reproducible/clinical | Strong nonlinear batch the reference never saw |
| Azimuth / Seurat anchor transfer | CCA/PCA anchors; supervised PCA projection | Yes (precomputed ref) | `prediction.score.max`, `mapping.score` | Multimodal refs (CITE-seq/WNN), curated tissue atlases, R shop | Filtered-anchor pathology when query is very divergent |
| scPoli | Conditional VAE + learnable sample embeddings + cell-type prototypes | Yes (built on scArches) | Prototype distance + uncertainty | Want sample-level (patient) embeddings too, many small batches | Few samples (condition embedding underdetermined) |
| popV | Ensemble of methods + ontology-aware voting | Mixed (wraps several) | Cross-method disagreement = uncertainty | High-stakes atlas annotation; distrust any single method | Compute-heavy; consensus can be confidently wrong if all share reference bias |
| CellTypist | Logistic regression (pre-trained models) | No embedding (ships models) | Low max-probability = ambiguous; no true OOD | Fast immune annotation, no integration needed | Treated as a mapper (no shared embedding, no batch handling) |
| treeArches / scHPL | scArches + hierarchical classifier with rejection | Yes | Explicit rejection -> "unseen" node | Expect novel subtypes, want hierarchy-aware fallback | Mis-specified hierarchy propagates error down branches |
| Foundation models (scGPT, Geneformer) | Transformer pretrained on 10s of millions of cells | Checkpoint; fine-tune needs labels | None intrinsic; OOD poorly characterized | *Fine-tuned* on the target task; cross-modality/species; data-scarce | Zero-shot: underperform scVI/Harmony/HVG-PCA (Kedzierska 2025) |

Cross-cutting: linear methods (Symphony, Azimuth-sPCA) give a *fixed, reproducible* reference embedding (the query never perturbs the reference); VAE surgery *fine-tunes* and can drift. Reproducibility-critical or clinical pipelines lean linear; maximal batch-effect flexibility leans VAE.

## Decision Tree by Scenario

| Scenario | Recommended approach | Why |
|----------|---------------------|-----|
| Same tissue as a published scVI/scANVI atlas; want one embedding + labels | scArches surgery onto the scANVI model; then gate labels on kNN uncertainty | Built for this; reuses the learned manifold; only query fine-tuned |
| Seconds-scale, deterministic, CPU-only, must re-run identically (clinical) | Symphony (or Azimuth if curated) | Fixed reference embedding, no fine-tuning drift, built-in Mahalanobis OOD |
| Query likely contains cell types/states NOT in the reference (disease, new niche) | treeArches/scHPL rejection, or scArches + explicit OOD distance | Default kNN/softmax confidently mislabels novel cells |
| High-stakes annotation; distrust any single method | popV ensemble | Cross-method disagreement is a more honest uncertainty than any softmax |
| Want patient/sample-level structure too | scPoli | Only method learning sample (condition) embeddings jointly with cell prototypes |
| Just need fast immune labels, no integration | CellTypist (`Immune_All_Low`, `majority_voting=True`) | Calibrated classifier, no embedding needed; QC the query first |
| Multimodal reference (CITE-seq/ATAC) | Azimuth/Seurat WNN or totalVI+scArches | Anchor framework natively weights modalities |
| Considering scGPT/Geneformer | Only if fine-tuning with labels, or cross-modality/species, or data too scarce | Zero-shot foundation embeddings are not a justified default for same-tissue transfer |
| De novo clustering, no reference, or manual marker annotation | -> single-cell/markers-annotation, single-cell/clustering | Out of scope here |

## scArches Surgery: Embedding (scVI)

**Goal:** Project query cells into a pre-trained reference latent space without retraining on combined data.

**Approach:** Align query genes to the reference exactly, load into the frozen reference model, and fine-tune only query-specific parameters with zero weight decay so the shared manifold does not drift.

```python
import scvi
import scanpy as sc

ref_model = scvi.model.SCVI.load('reference_model/')          # saved with save_anndata=True (or minified)

adata_query = sc.read_h5ad('query.h5ad')
# Align genes to the reference EXACTLY: zero-pad missing, reorder. Mandatory and silent if skipped.
scvi.model.SCVI.prepare_query_anndata(adata_query, 'reference_model/')

query_model = scvi.model.SCVI.load_query_data(adata_query, 'reference_model/')
# weight_decay=0.0 + frozen reference weights make surgery a query-only fine-tune;
# non-zero decay drifts the shared latent and breaks cross-query comparability.
query_model.train(max_epochs=200, plan_kwargs={'weight_decay': 0.0}, check_val_every_n_epoch=10)
adata_query.obsm['X_scVI'] = query_model.get_latent_representation()
```

## scANVI Label Transfer

**Goal:** Transfer reference cell-type labels to an unlabeled query.

**Approach:** Build a semi-supervised scANVI head on the reference, map the query by surgery, then read hard labels and per-class probabilities -- treating the probability as "which label," not "does it belong."

```python
# Reference side (once): scANVI from a trained scVI model. unlabeled_category is REQUIRED.
ref_scanvi = scvi.model.SCANVI.from_scvi_model(ref_vae, unlabeled_category='Unknown', labels_key='cell_type')
ref_scanvi.train(max_epochs=20, n_samples_per_label=100)
ref_scanvi.save('ref_scanvi/', save_anndata=True)

# Query side (surgery):
scvi.model.SCANVI.prepare_query_anndata(adata_query, 'ref_scanvi/')
query_scanvi = scvi.model.SCANVI.load_query_data(adata_query, 'ref_scanvi/')
query_scanvi.train(max_epochs=100, plan_kwargs={'weight_decay': 0.0})

adata_query.obs['predicted_label'] = query_scanvi.predict()          # hard labels
adata_query.obsm['X_scANVI'] = query_scanvi.get_latent_representation()
proba = query_scanvi.predict(soft=True)                             # per-class probabilities (which label)
```

## Out-of-Distribution Gating (the step that makes labels trustworthy)

**Goal:** Decide *whether* each query cell belongs to the reference, separately from which label it would get.

**Approach:** Compute a distance/entropy signal on the shared latent. The canonical scArches/HLCA approach is a weighted-kNN label-transfer uncertainty (neighbor disagreement in the reference latent), thresholded at 0.2 to set cells to "Unknown." A portable kNN-entropy version is shown; the softmax `proba` is NOT this signal.

```python
import numpy as np
from sklearn.neighbors import KNeighborsClassifier

ref_latent = ref_scanvi.get_latent_representation()                  # reference cells in latent
knn = KNeighborsClassifier(n_neighbors=15, weights='distance').fit(ref_latent, adata_ref.obs['cell_type'])
query_latent = adata_query.obsm['X_scANVI']

neighbor_proba = knn.predict_proba(query_latent)                    # weighted neighbor label distribution
# Uncertainty = 1 - max neighbor agreement. HLCA sets cells above 0.2 to 'Unknown'.
uncertainty = 1.0 - neighbor_proba.max(axis=1)
adata_query.obs['transfer_uncertainty'] = uncertainty
adata_query.obs.loc[uncertainty > 0.2, 'predicted_label'] = 'Unknown'   # gate, do not trust ungated labels
print(f'Flagged Unknown: {(uncertainty > 0.2).mean():.1%}')
```

## Per-Method Failure Modes

### scANVI / kNN -- forcing the query onto reference labels
- **Trigger:** Query contains a population absent from the reference (novel type, disease state, perturbed program).
- **Mechanism:** The classifier/kNN assigns every query cell to its nearest reference label; there is no "none of the above" unless added.
- **Symptom:** A coherent novel cluster split across 2-3 reference labels, each with *high* probability; the UMAP looks "well integrated."
- **Fix:** Always gate on transfer uncertainty / OOD distance (above). Inspect query-only clusters for marker genes independent of transferred labels.

### Softmax overconfidence vs label-transfer uncertainty conflated
- **Trigger:** Reporting "confidence" as the `predict(soft=True)` max.
- **Mechanism:** The softmax measures *which* reference label conditional on belonging; it is normalized away from distance and cannot say "far from everything." A cell can be 0.99 "T cell" and be a hepatocyte.
- **Symptom:** Pipeline filters on softmax >= 0.5 and still passes OOD cells.
- **Fix:** Threshold the weighted-kNN uncertainty (HLCA 0.2) or a Mahalanobis/ensemble OOD signal for the "Unknown" decision; use the softmax only to disambiguate among in-distribution labels.

### scANVI -- semi-supervised label leakage / latent carving
- **Trigger:** scANVI reference where labels strongly drive latent geometry; trajectory or novel-state query.
- **Mechanism:** The classifier head back-propagates label structure into the latent, carving it to separate reference types; query cells are pulled toward that structure even when their biology lies between/outside it.
- **Symptom:** A query continuum (differentiation trajectory) collapses onto discrete reference clusters; intermediate states vanish.
- **Fix:** For trajectory/novel-state queries prefer *unsupervised* scVI surgery, annotate the query independently, and cross-check against the scVI latent.

### Reference composition / missing-biology bias
- **Trigger:** Reference from healthy/limited donors; query from disease, different age, ancestry, or tissue region.
- **Mechanism:** The reference manifold is the entire hypothesis space; off-manifold cells are projected onto the nearest on-manifold point and rare reference populations act as attractors.
- **Symptom:** Disease-specific states labeled as the closest healthy type; ancestry/age effects read as "batch."
- **Fix:** Audit reference composition before mapping; prefer references covering the query's expected biology; treat mapping as hypothesis generation and validate query findings de novo; consider extending the reference (treeArches).

### Query QC artifacts laundered into confident labels
- **Trigger:** Query not QC'd to the reference's standard (empty droplets, ambient RNA, doublets, high-MT).
- **Mechanism:** A doublet sits between two reference types and maps to a spurious "intermediate"; ambient RNA shifts profiles toward the dominant type.
- **Symptom:** Artifactual "transitional" populations; doublet clusters labeled as rare real types.
- **Fix:** Run full query QC *before* mapping (single-cell/doublet-detection, ambient correction, MT/count filters matched to the reference). Mapping does not clean data.

### Feature-space / gene-set mismatch
- **Trigger:** Query missing reference HVGs; different gene annotation/version; `prepare_query_anndata` skipped.
- **Mechanism:** The encoder expects the exact reference gene order; missing genes are zero-padded and reordering silently corrupts the input.
- **Symptom:** Garbage latent, everything maps to one blob, or a silent accuracy cliff (no error raised).
- **Fix:** Always `prepare_query_anndata(query, reference_model)`; verify the shared-gene fraction; too few shared HVGs is a hard stop.

### Good integration metrics, wrong labels
- **Trigger:** Judging mapping by scIB integration scores alone.
- **Mechanism:** Integration metrics reward query/reference mixing; mixing OOD cells into the wrong neighborhood *raises* the batch-removal score, and bio-conservation uses reference labels (circular for the query).
- **Symptom:** Top scIB total score with biologically wrong annotation.
- **Fix:** Integration metrics validate the *embedding*, not labels. Evaluate transfer on held-out labeled query cells (per-type F1, especially rare types), OOD detection on spiked-in unseen types, and marker-gene sanity checks.

### Zero-shot foundation-model embedding as a mapper
- **Trigger:** Using scGPT/Geneformer zero-shot embeddings for clustering/transfer expecting SOTA.
- **Mechanism:** The masked-gene pretraining objective does not guarantee a label- or batch-aware embedding; zero-shot embeddings are not batch-corrected.
- **Symptom:** Worse AvgBio/integration than scVI or even HVG-PCA + Harmony.
- **Fix:** Fine-tune with task labels, or use an established mapper; reserve foundation models for cross-modality/species/data-scarce cases (Kedzierska 2025).

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| scANVI label confident but kNN uncertainty high | OOD cell forced onto nearest label | Trust the uncertainty; set Unknown and inspect markers |
| Symphony Mahalanobis flags OOD but scANVI does not | scANVI latent carved to absorb the cell | Prefer the distance-based flag; novel biology likely |
| popV members disagree | Genuine ambiguity or granularity mismatch | Route to manual review; report the disagreement, do not force a leaf |
| High scIB score, poor per-type F1 on held-out labels | Embedding mixes well but labels wrong | Believe the F1; integration score is not a label metric |

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Transfer uncertainty > 0.2 -> "Unknown" | Sikkema 2023 (HLCA) | Weighted-kNN neighbor-disagreement cutoff bounding false labels; recalibrate per reference |
| Surgery `weight_decay=0.0`, ~100-200 epochs | scvi-tools scArches tutorial | Frozen reference weights + no decay keep the shared latent fixed |
| scIB total = 0.6*bio + 0.4*batch | Luecken 2022 | Benchmark weighting; scores the embedding, NOT query labels |
| CellTypist input = log1p of CP10k | CellTypist docs | Wrong normalization silently degrades accuracy |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Everything maps to one blob | `prepare_query_anndata` skipped; gene mismatch | Run it before `load_query_data`; check shared-gene fraction |
| OOD cells pass a 0.5 softmax filter | Thresholding the wrong quantity | Gate on weighted-kNN uncertainty / OOD distance, not softmax |
| Reference cells move between runs | Non-zero `weight_decay` or over-training in surgery | Set `weight_decay=0.0`, keep `freeze_*` defaults, modest epochs |
| `load_query_data` errors on labels | scANVI needs `unlabeled_category` | Pass it to `from_scvi_model`; query labels filled with that category |
| CellTypist labels look random | Raw or wrongly normalized counts | Feed log1p CP10k input |

## References

- Lopez R, Regier J, Cole MB, Jordan MI, Yosef N. 2018. Deep generative modeling for single-cell transcriptomics. *Nat Methods* 15:1053-1058.
- Xu C, Lopez R, Mehlman E, Regier J, Jordan MI, Yosef N. 2021. Probabilistic harmonization and annotation of single-cell transcriptomics data with deep generative models. *Mol Syst Biol* 17:e9620.
- Lotfollahi M, Naghipourfar M, Luecken MD, et al. 2022. Mapping single-cell data to reference atlases by transfer learning. *Nat Biotechnol* 40:121-130.
- Kang JB, Nathan A, Weinand K, et al. 2021. Efficient and precise single-cell reference atlas mapping with Symphony. *Nat Commun* 12:5890.
- Hao Y, Hao S, Andersen-Nissen E, et al. 2021. Integrated analysis of multimodal single-cell data. *Cell* 184:3573-3587.
- De Donno C, Hediyeh-Zadeh S, Moinfar AA, et al. 2023. Population-level integration of single-cell datasets enables multi-scale analysis across samples. *Nat Methods* 20:1683-1692.
- Ergen C, Xing G, Xin C, et al. 2024. Consensus prediction of cell type labels in single-cell data with popV. *Nat Genet* 56:2731-2738.
- Dominguez Conde C, Xu C, Jarvis LB, et al. 2022. Cross-tissue immune cell analysis reveals tissue-specific features in humans. *Science* 376:eabl5197.
- Michielsen L, Lotfollahi M, Strobl D, et al. 2023. Single-cell reference mapping to construct and extend cell-type hierarchies. *NAR Genom Bioinform* 5:lqad070.
- Luecken MD, Buttner M, Chaichoompu K, et al. 2022. Benchmarking atlas-level data integration in single-cell genomics. *Nat Methods* 19:41-50.
- Sikkema L, Ramirez-Suastegui C, Strobl DC, et al. 2023. An integrated cell atlas of the lung in health and disease. *Nat Med* 29:1563-1577.
- Kedzierska KZ, Crawford L, Amini AP, Lu AX. 2025. Zero-shot evaluation reveals limitations of single-cell foundation models. *Genome Biol* 26:101.

## Related Skills

- single-cell/preprocessing - QC, normalization, and HVG selection the query needs before mapping
- single-cell/markers-annotation - Manual marker-based cluster annotation when there is no reference
- single-cell/batch-integration - Integrating datasets without a labeled reference
- single-cell/doublet-detection - Removing doublets that map to spurious intermediates
- differential-expression/de-results - Pseudobulk validation of mapping-derived populations
