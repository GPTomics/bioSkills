---
name: bio-single-cell-multimodal-integration
description: Integrate multimodal single-cell data (CITE-seq RNA+protein, 10x Multiome RNA+ATAC, unpaired/diagonal RNA+ATAC) and choose the right joint method. Use when classifying an integration task by anchor structure (paired vs unpaired), denoising CITE-seq ADT background before joint embedding, picking between WNN, totalVI, MultiVI, MOFA+, GLUE, or Seurat v5 bridge integration, or diagnosing why a modality dominates a joint clustering.
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

# Multimodal Integration

**"Jointly analyze my CITE-seq / Multiome / unpaired multi-omic data"** -> Classify the task by anchor structure, denoise each modality in its native pipeline, then build one joint representation.
- R: `Seurat::FindMultiModalNeighbors()` (WNN), `Signac` (ATAC LSI), `dsb::DSBNormalizeProtein()` (ADT denoising), `PrepareBridgeReference()` (v5 bridge)
- Python: `muon`/`mudata` (MuData container), `scvi.model.TOTALVI` / `MULTIVI`, `MOFA2`/`muon.tl.mofa`, `scglue` (diagonal)

## Governing Principle

Classify the integration task by its anchor structure FIRST, because the anchor decides which algorithm class is even applicable (Argelaguet 2021).
- Horizontal: same modality, different cells; anchor = shared features (batch correction, not this skill).
- Vertical (paired, same cell): multiple modalities measured in the same cells; anchor = shared cells. CITE-seq, 10x Multiome.
- Diagonal (unpaired): different modalities in different cells, no shared cells and no shared features; correspondence is inferred from prior knowledge. Independent scRNA + scATAC.
- Mosaic: partially observed grid of (modalities x batches); some blocks present, some missing.

Paired vs unpaired is the master fork: paired correspondence is known a priori (WNN, totalVI, MultiVI, MOFA+), unpaired/diagonal correspondence must be inferred (GLUE, Seurat v5 bridge), mosaic mixes both (MultiVI, StabMap). Two separately-paired datasets that share only one modality (for example a 10x Multiome and a CITE-seq experiment sharing only RNA) are a mosaic problem: anchor on the shared RNA and impute or bridge the modality-specific blocks with StabMap, MultiVI, or Seurat v5 bridge integration rather than forcing a single WNN.

CITE-seq ADT background is a three-part mixture, not one "ambient" term: (1) ambient antibody captured in every droplet including empties, (2) cell-intrinsic non-specific binding (Fc receptors, sticky dying cells) that does NOT appear in empties, (3) spillover/index hopping between barcodes. Denoise ADT (DSB or totalVI's built-in background mixture) BEFORE any joint embedding; raw or CLR-only ADT carries this background into the joint graph.

WNN can be dominated by the noisier modality: weights reward local neighbor predictability, and a handful of high-variance or saturating ADT features can manufacture self-consistent neighborhoods and get up-weighted despite carrying less biology. Report the per-cell weight distribution and check whether clustering survives down-weighting the suspect modality.

Imputed modalities are inferences, not measurements: MultiVI/StabMap/Cobolt impute the missing modality for unpaired cells, and gene-activity scores from ATAC approximate RNA. Differential expression or marker calls on imputed values are model-dependent and must be flagged as such.

## Classify the Task: Anchor Structure -> Method Class

| Anchor structure | What is shared | Example assay | Method class |
|---|---|---|---|
| Vertical / paired | Same cells | CITE-seq, 10x Multiome | WNN, totalVI, MultiVI(paired), MOFA+, mojitoo |
| Diagonal / unpaired | Nothing (prior graph) | Independent scRNA + scATAC | GLUE, Seurat v5 bridge, LIGER |
| Mosaic | Some modalities only | Batch A RNA+ATAC, batch B RNA | MultiVI, StabMap, Cobolt, totalVI(partial) |

When methods compete, verify the current best-practice default against the installed tool docs before committing; the field moves and defaults drift across minor versions.

## Method Decision Table (Paired CITE-seq / Multiome)

| Method | Model / assumption | Use when | Fails when |
|---|---|---|---|
| WNN (Seurat) | Per-cell, per-modality weights from cross-modality neighbor prediction; one weighted graph | Fast joint embedding/clustering of one well-normalized paired dataset | Protein background not removed upstream; noisy/saturating modality dominates; not for unpaired/mosaic |
| totalVI (scvi-tools) | Conditional VAE; RNA NB/ZINB, each protein a 2-component NB mixture (background+foreground) | Need denoised protein, principled DE, batch integration, merging different antibody panels | Tiny datasets (VAE overfits); no GPU and very large data; protein-specific background structure not captured by one per-cell factor |
| MultiVI (scvi-tools) | Single joint VAE over RNA+ATAC(+protein); mosaic-capable, imputes missing modality | Paired+unpaired RNA/ATAC mixed (mosaic); want generative DE/DA | "batch" key is the modality indicator, not sequencing batch; imputed modalities treated as measured |
| MOFA+ (MOFA2) | Linear Bayesian group factor analysis; sparse factors, per-modality variance explained | Interpreting shared vs modality-specific axes of variation (exploratory/explanatory) | Used for clustering/denoising; likelihood mismatched to data; expecting batch correction within a view |
| mojitoo | CCA across precomputed per-modality reductions; fast, parameter-free | Quick paired joint reduction from existing PCA/LSI slots | No knob to down-weight a noisy modality; bounded by input reductions; paired only |

## Method Decision Table (Unpaired / Diagonal / Mosaic)

| Method | Model / assumption | Use when | Fails when |
|---|---|---|---|
| GLUE (scglue) | Per-modality VAEs + prior feature graph (peak-near-gene); adversarial cell alignment | Unpaired diagonal scRNA + scATAC; want regulatory inference as a byproduct | Genome-build/coordinate mismatch yields an empty guidance graph and garbage alignment; adversarial over-mixing of distinct states |
| Seurat v5 bridge | Multiome bridge dataset = dictionary linking query modality to reference modality | Mapping a query (scATAC) onto a reference built in another modality (scRNA) | Poor/batch-mismatched bridge propagates error; rare query-only populations mislabeled |
| StabMap | Mosaic topology from shared features; project all cells via shortest paths | Mosaic with informative unshared features that cannot be dropped | Unshared-feature chaining compounds error per hop |
| Cobolt / scMoMaT | Generative shared latent over joint + single-modality datasets | Mosaic where a generative latent is preferred over feature chaining | DE/marker calls made on imputed values |

## ADT Normalization: CLR vs DSB

| Method | What it does | Use when | Fails when |
|---|---|---|---|
| CLR (centered log-ratio) | Rescales compositionally; Seurat `NormalizeData(method="CLR", margin=2)` | Quick, no empty droplets available; small panels | Does NOT remove background; geometric-mean denominator distorted by saturating high-abundance ADTs |
| DSB | Ambient correction from empty droplets + per-cell technical denoising via 2-component mixture + isotype controls | Raw/unfiltered matrix available (needs empty droplets); want background removed before embedding | No empty droplets retained; protein-specific non-specific binding (one per-cell factor under/over-corrects); no clearly bimodal proteins |

Seurat's CLR margin is genuinely ambiguous across versions (margin=2 = per-feature is the WNN-tutorial recommendation for large panels); verify with `?NormalizeData` on the installed version.

## CITE-seq: Denoise ADT, Then Joint Embed (Seurat)

**Goal:** Remove ADT background with DSB before WNN, because WNN does not denoise protein.

**Approach:** Estimate ambient from empty droplets and per-cell technical noise from a mixture plus isotype controls, then feed denoised ADT into the standard PCA -> WNN flow.

```r
library(dsb)
library(Seurat)

raw <- Read10X('raw_feature_bc_matrix/')           # unfiltered: contains empty droplets
cells <- Read10X('filtered_feature_bc_matrix/')    # called cells

adt_cells <- as.matrix(cells[['Antibody Capture']])
adt_empty <- as.matrix(raw[['Antibody Capture']][, setdiff(colnames(raw[['Antibody Capture']]), colnames(adt_cells))])

# isotype.control.name.vec must name the ACTUAL isotype rows (often IgG1/IgG2a/Mouse-IgG2b-Ctrl); the regex below misses those
# When isotypes are absent or not matched, set use.isotype.control = FALSE (keep denoise.counts = TRUE) and pass real names explicitly
adt_dsb <- DSBNormalizeProtein(
    cell_protein_matrix = adt_cells,
    empty_drop_matrix = adt_empty,
    denoise.counts = TRUE,
    use.isotype.control = TRUE,
    isotype.control.name.vec = grep('[Ii]sotype|IgG', rownames(adt_cells), value = TRUE)
)
```

## CITE-seq: WNN Joint Clustering (Seurat)

**Goal:** Build one weighted-NN graph from denoised RNA and ADT and cluster on it.

**Approach:** Reduce each modality independently (PCA on RNA, PCA on the small ADT panel), then learn per-cell modality weights and cluster/embed on the joint graph.

```r
obj[['ADT']] <- CreateAssay5Object(data = adt_dsb)        # DSB output is already normalized data
DefaultAssay(obj) <- 'RNA'
obj <- NormalizeData(obj) |> FindVariableFeatures() |> ScaleData() |> RunPCA(reduction.name = 'pca')

DefaultAssay(obj) <- 'ADT'
VariableFeatures(obj) <- rownames(obj[['ADT']])
obj <- ScaleData(obj) |> RunPCA(reduction.name = 'apca', npcs = min(18, nrow(obj[['ADT']]) - 1))

# dims.list matched to informative dims; small ADT panels saturate by ~1:18
obj <- FindMultiModalNeighbors(obj, reduction.list = list('pca', 'apca'), dims.list = list(1:30, 1:18))
obj <- FindClusters(obj, graph.name = 'wsnn', algorithm = 3)   # algorithm 3 = SLM (the tutorial choice), NOT Leiden
obj <- RunUMAP(obj, nn.name = 'weighted.nn', reduction.name = 'wnn.umap')

# Inspect the per-cell weight distribution; a single dominant modality is a red flag
VlnPlot(obj, features = 'RNA.weight', group.by = 'seurat_clusters')
```

## CITE-seq: totalVI (Python, denoise + DE in one model)

**Goal:** Jointly model RNA + protein with explicit protein background, yielding a denoised latent space and foreground probabilities.

**Approach:** Register a MuData object, train the conditional VAE, then read the latent representation and per-protein foreground probability.

```python
import scvi
import mudata as md

# mdata holds .mod['rna'] (raw counts) and .mod['prot'] (raw ADT counts)
scvi.model.TOTALVI.setup_mudata(
    mdata, rna_layer='counts', protein_layer=None,
    modalities={'rna_layer': 'rna', 'protein_layer': 'prot'}
)
model = scvi.model.TOTALVI(mdata)
model.train()

mdata.obsm['X_totalVI'] = model.get_latent_representation()
fg = model.get_protein_foreground_probability()        # 1 - background mixing weight per protein per cell
denoised_rna, denoised_prot = model.get_normalized_expression()
```

## Multiome (RNA + ATAC, same cell): Native Pipelines, Then Join

**Goal:** Process each modality in its own statistics before joining, because RNA and ATAC have incompatible distributions.

**Approach:** PCA on RNA, TF-IDF + LSI on ATAC (drop depth-correlated components), then WNN. See scatac-analysis for ATAC QC and the binarization/depth-component caveats.

```r
library(Signac)
DefaultAssay(obj) <- 'RNA'
obj <- NormalizeData(obj) |> FindVariableFeatures() |> ScaleData() |> RunPCA()

DefaultAssay(obj) <- 'ATAC'
obj <- RunTFIDF(obj) |> FindTopFeatures(min.cutoff = 'q0') |> RunSVD()
DepthCor(obj)                                          # diagnose which LSI components track depth

# dims = 2:30 drops LSI_1 ONLY if DepthCor confirms it tracks depth (usually true, not guaranteed)
obj <- FindMultiModalNeighbors(obj, reduction.list = list('pca', 'lsi'), dims.list = list(1:30, 2:30))
obj <- RunUMAP(obj, nn.name = 'weighted.nn', reduction.name = 'wnn.umap')
obj <- FindClusters(obj, graph.name = 'wsnn', algorithm = 3)
```

Merging multiome datasets requires a common peak set: re-quantify all cells against unified peaks, or peak-boundary differences manufacture spurious batch structure. The ATAC gene-activity matrix is an approximation, not measured RNA; do not conflate it with the RNA modality.

## MOFA+ (interpretable shared/specific factors)

**Goal:** Decompose modalities into shared latent factors with per-modality variance explained.

**Approach:** Build a MOFA object from per-modality matrices, set likelihoods to match each data type, run, then interpret factor loadings.

```python
import muon as mu

# likelihoods must match data: gaussian for scaled RNA, bernoulli for binarized ATAC, poisson for counts
mu.tl.mofa(mdata, n_factors=15, outfile='mofa_model.hdf5')   # writes mdata.obsm['X_mofa']
```

## Unpaired / Diagonal: GLUE (Python)

**Goal:** Align independent scRNA and scATAC with no shared cells via a prior feature graph.

**Approach:** Configure each dataset with a count-appropriate probabilistic model, build a gene-anchored guidance graph, fit GLUE, then read aligned embeddings.

```python
import scglue

scglue.models.configure_dataset(rna, 'NB', use_highly_variable=True, use_rep='X_pca')     # NB needs RAW counts
scglue.models.configure_dataset(atac, 'ZINB', use_highly_variable=True, use_rep='X_lsi')
graph = scglue.genomics.rna_anchored_guidance_graph(rna, atac)     # peak-near-gene prior; coords must share genome build
glue = scglue.models.fit_SCGLUE({'rna': rna, 'atac': atac}, graph)
rna.obsm['X_glue'] = glue.encode_data('rna', rna)
atac.obsm['X_glue'] = glue.encode_data('atac', atac)
```

Verify cell-type structure is preserved (not just modality overlap); adversarial alignment can over-mix distinct populations.

## MuData Housekeeping

After per-modality QC, modalities hold different cell sets; `muon.pp.intersect_obs(mdata)` before any paired analysis. Editing a modality-local `mdata.mod['rna'].obs` needs `mdata.update()` to propagate to the global `mdata.obs`. R round-trips (MuDataSeurat, zellkonverter) are lossy; plan to stay in one ecosystem.

## Common Errors

| Symptom | Cause | Fix |
|---|---|---|
| WNN clustering driven entirely by ADT | A few saturating high-variance proteins dominate the neighbor graph | Report per-cell weight distribution; down-weight or denoise ADT (DSB); re-check clustering stability |
| "Background" smear in every ADT cluster | Ran WNN/CLR without empty-droplet denoising | Run DSB (needs raw/unfiltered matrix) or totalVI before joint embedding |
| DSB errors / nonsense output | Passed a filtered cell matrix only (no empty droplets) | Supply `empty_drop_matrix` from the raw/unfiltered matrix |
| Spurious batch structure after merging multiome | Per-dataset peak sets, not a unified set | Re-quantify all cells against one common peak set |
| GLUE produces a blob / no alignment | Guidance graph near-empty from genome-build/coordinate mismatch | Align RNA gene coords and ATAC peaks to the same build before building the graph |
| RNA and protein disagree for a marker | Often real post-transcriptional biology (stability, trafficking, lag), not an artifact | Do not "correct away"; treat single-gene discordance as informative |
| MultiVI batch effects persist | The `batch_key` was set to the modality indicator, not sequencing batch | Add a separate covariate for the real batch |
| DE on a modality looks too clean | Computed on imputed/gene-activity values, not measurements | Flag imputed-modality DE as model-dependent; validate against a measured modality |

## Related Skills

- single-cell/scatac-analysis - ATAC QC, TF-IDF/LSI, gene-activity caveats for the Multiome ATAC half
- single-cell/preprocessing - per-modality RNA QC and normalization before integration
- single-cell/clustering - clustering and UMAP on the joint graph
- single-cell/batch-integration - horizontal (same-modality, cross-sample) correction
- single-cell/markers-annotation - marker-based interpretation of joint clusters
- atac-seq/motif-deviation - chromVAR TF activity on the Multiome ATAC modality
- pathway-analysis/go-enrichment - functional interpretation of modality-specific factors

## References

Argelaguet R, Cuomo ASE, Stegle O, Marioni JC. Computational principles and challenges in single-cell data integration. Nat Biotechnol 39(10):1202-1215 (2021).
Stoeckius M, Hafemeister C, Stephenson W, et al. Simultaneous epitope and transcriptome measurement in single cells (CITE-seq). Nat Methods 14:865-868 (2017).
Mulè MP, Martins AJ, Tsang JS. Normalizing and denoising protein expression data from droplet-based single-cell profiling (DSB). Nat Commun 13:2099 (2022).
Hao Y, Hao S, Andersen-Nissen E, et al. Integrated analysis of multimodal single-cell data (WNN). Cell 184(13):3573-3587 (2021).
Gayoso A, Steier Z, Lopez R, et al. Joint probabilistic modeling of single-cell multi-omic data with totalVI. Nat Methods 18:272-282 (2021).
Ashuach T, Gabitto MI, Koodli RV, et al. MultiVI: deep generative model for the integration of multimodal data. Nat Methods 20(8):1222-1231 (2023).
Argelaguet R, Arnol D, Bredikhin D, et al. MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data. Genome Biol 21:111 (2020).
Cao Z-J, Gao G. Multi-omics single-cell data integration and regulatory inference with graph-linked unified embedding (GLUE). Nat Biotechnol 40(10):1458-1466 (2022).
Hao Y, Stuart T, Kowalski MH, et al. Dictionary learning for integrative, multimodal and scalable single-cell analysis (Seurat v5 bridge). Nat Biotechnol 42:293-304 (2024).
Bredikhin D, Kats I, Stegle O. MUON: multimodal omics analysis framework. Genome Biol 23:42 (2022).
Ghazanfar S, Guibentif C, Marioni JC. Stabilized mosaic single-cell data integration using unshared features (StabMap). Nat Biotechnol 42(2):284-292 (2024).
Yin Y, et al. Characterization and decontamination of background noise in droplet-based single-cell protein expression data with DecontPro. Nucleic Acids Res 52(1):e4 (2024).
