---
name: bio-single-cell-cell-communication
description: Infers ligand-receptor cell-cell communication from scRNA-seq with a consensus-first workflow (LIANA), plus CellPhoneDB specificity tests, CellChat pathway probabilities, and NicheNet downstream ligand-activity. Use when ranking ligand-receptor interactions between cell types, comparing communication across conditions, asking which ligand drives a receiver response, or deciding which CCC method and resource to trust.
tool_type: mixed
primary_tool: LIANA
---

## Version Compatibility

Reference examples tested with: liana 1.2+, CellChat 2.1+, cellphonedb 5.0+, nichenetr 2.1+, scanpy 1.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Cell-Cell Communication Analysis

**"Find which cell types signal to each other"** -> Score ligand-receptor pairs from co-expression in sender and receiver populations, rank them, and assess specificity or downstream effect.
- Python: `liana.mt.rank_aggregate()` (consensus default), `cellphonedb` (permutation specificity)
- R: `CellChat::computeCommunProb()` (pathway probability), `nichenetr::predict_ligand_activities()` (downstream mechanism)

## Governing Principle

Every ligand-receptor output is a co-expression PROXY, not proof of signaling. Co-expression is neither necessary nor sufficient: receptor mRNA is not a responsive surface protein (desensitization, internalization, decoy receptors, missing co-receptors), a ligand must be secreted or proteolytically cleaved and activated to act, and spatial proximity is unobserved in dissociated data so the two cell types may never have been adjacent. Methods are DISCORDANT because they estimate DIFFERENT quantities, not because of noise: CellPhoneDB tests expression SPECIFICITY (label permutation), CellChat tests a mass-action PROBABILITY (magnitude + Hill saturation), NATMI/Connectome score MAGNITUDE, and there is no monotone mapping between these rankings, so top-N lists genuinely differ on identical input. The choice of RESOURCE (the L-R database) can move results as much as or more than the choice of method (Dimitrov 2022). Default to a CONSENSUS (LIANA `rank_aggregate`), run a resource-sensitivity check, treat every interaction as a hypothesis, and validate orthogonally (downstream TF/pathway activity, receptor protein by CITE-seq, spatial co-localization, or perturbation). NicheNet asks a distinct, better-grounded question (which ligand explains the receiver's DE response) but depends on a clean receiver gene set and a static, cell-type-agnostic prior network.

## Method Decision Table

| Method | Tests what / null | Use when | Fails when |
|--------|-------------------|----------|------------|
| LIANA `rank_aggregate` | Consensus rank over many scoring functions; reports magnitude AND specificity ranks | Robust default; hedge against discordance; vary resource to test sensitivity | Treated as ground truth; top-N read as stable (tail of ranks is flat, membership is unstable to subsampling/re-clustering) |
| CellPhoneDB v5 | Expression SPECIFICITY; permutes cluster labels, asks if mean L+R expression exceeds random labelling | Permutation p-values, rigorous multi-subunit complexes (limiting subunit), human, spatial microenvironments / CellSign TF add-on | Mouse data (human-only DB, ortholog mapping errors); abundance drives the null so dominant clusters over-call; magnitude ignored |
| CellChat v2 | Communication PROBABILITY; law-of-mass-action + Hill saturation, cofactor terms, trimean expression | Pathway-level summaries, sender/receiver/mediator roles, cofactor modeling, cross-condition comparison, fewer high-confidence calls | Sparse/lowly expressed genes dropped by conservative trimean; interaction COUNTS compared across datasets without normalization |
| NATMI / Connectome | MAGNITUDE; expression product (NATMI adds a specificity edge weight) | A simple, fast magnitude score; component of LIANA consensus | Used alone as "communication" - pure magnitude rewards ubiquitous high genes |
| NicheNet | Downstream LIGAND-ACTIVITY; ranks ligands by AUPR between predicted regulatory targets and the receiver's observed DE genes | The question is mechanism: which ligand best explains THIS receiver response | Receiver gene set is noisy/batch-confounded; prior network is static and cell-type-agnostic so context-specific wiring is missed; not a de-novo who-talks-to-whom tool |

Methods evolve; before committing, verify current best practice and the default resource against the installed package docs (LIANA NEWS, CellChatDB version, cellphonedb-data release).

For communication PROGRAMS varying across many samples/conditions/time, decompose with Tensor-cell2cell (commonly run as LIANA -> Tensor-cell2cell) rather than comparing raw counts; for comparison at single-cell resolution without cluster averaging, use Scriabin, which recovers edges lost to agglomeration.

## Spatial-Aware Methods (proximity != interaction)

In dissociated scRNA-seq, proximity is UNKNOWN - only spatial methods constrain by physical distance, and even then they demonstrate spatially-coherent co-expression under modeling assumptions, not binding.

| Method | Approach | Caveat |
|--------|----------|--------|
| Squidpy `sq.gr.ligrec` | CellPhoneDB-style permutation on spatial coordinates | Visium spots hold multiple cells -> "co-expression" can be two cells in one spot; deconvolve first |
| COMMOT | Collective optimal transport over a distance cost with a hard diffusion radius | Radius and cost kernel are unvalidated hyperparameters; run a sensitivity analysis over the radius |
| CellChat v2 spatial | Mass-action probability constrained by spatial distance | Same trimean conservatism; distance scaling is a modeling choice |
| LIANA+ bivariate | Local bivariate (Moran's-style) L-R co-occurrence in space | Tests spatial co-distribution, confounded by shared niche regulation; not directed signaling |

## Confounds That Mimic Signaling

| Confound | How it manufactures a fake interaction | Mitigation |
|----------|----------------------------------------|------------|
| Ambient RNA | Soup of highly expressed secreted genes (hemoglobin, albumin, cytokines) leaks into every cluster, inflating the SECRETED-ligand half of pairs and creating "universal senders" | Decontaminate (SoupX/DecontX/CellBender) before CCC, especially for secreted ligands |
| Cell-type abundance | Larger clusters give a tighter permutation null and smaller p-values, so the dominant type becomes the hub of every network | Down-sample or check that interaction counts do not simply track cluster sizes |
| Sequencing depth | Depth differences across samples change detected genes and scores, conflating technical with biological signal | Run on integrated counts; normalize depth before cross-condition comparison |
| Dissociation stress | Enzymatic dissociation induces FOS, JUN, JUNB, EGR1, HSPA1A/B, DUSP1 - several are bona fide ligands, fabricating AP-1 / heat-shock "signaling" | Flag or regress the stress-gene module; consider cold-protease or snRNA-seq |

## Consensus Inference (Default)

**Goal:** Rank ligand-receptor pairs robustly without committing to one method's estimand.

**Approach:** Run LIANA's rank aggregation over many scoring functions on one input and one resource, then read BOTH the magnitude and specificity ranks (a pair can score high on one and low on the other). CCC needs >=2 cell types in `groupby`; a single group yields only autocrine self-edges, not intercellular signaling.

```python
import liana as li
import scanpy as sc

adata = sc.read_h5ad('adata_annotated.h5ad')

# expr_prop=0.1 drops pairs expressed in <10% of a cluster (sparse-noise floor)
# n_perms=1000 builds the specificity null; use_raw=False uses log-normalized .X
li.mt.rank_aggregate(adata, groupby='cell_type', resource_name='consensus',
                     expr_prop=0.1, use_raw=False, n_perms=1000, verbose=True)

res = adata.uns['liana_res']
# rank_aggregate yields magnitude_rank and specificity_rank (NOT a single 'liana_rank')
robust = res[(res['specificity_rank'] < 0.05) & (res['magnitude_rank'] < 0.05)]
```

## Resource-Sensitivity Check

**Goal:** Establish that a finding is not an artifact of one L-R database.

**Approach:** Hold the method fixed and re-run with a second resource; a pair that survives both resources is robust, one that flips is not.

```python
from liana.method import cellphonedb

for resource in ['consensus', 'cellphonedb', 'cellchatdb']:
    cellphonedb(adata, groupby='cell_type', resource_name=resource,
                expr_prop=0.1, use_raw=False, key_added=f'cpdb_{resource}', verbose=False)
# Compare top pairs across adata.uns['cpdb_consensus'] / 'cpdb_cellphonedb' / 'cpdb_cellchatdb'
```

## Specificity Test (CellPhoneDB v5)

**Goal:** Get permutation specificity p-values with rigorous multi-subunit complex handling (human).

**Approach:** Run the statistical method on log-normalized counts plus a cell-type meta table; the permutation null shuffles cluster labels, and complexes require all subunits via the limiting (minimum) subunit.

```python
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

# threshold=0.1: a gene must be expressed in >=10% of a cluster's cells to count
# iterations=1000: label-permutation null; pvalue=0.05 reports per-pair significance
results = cpdb_statistical_analysis_method.call(
    cpdb_file_path='cellphonedb.zip',          # cellphonedb-data v5 release
    meta_file_path='meta.tsv',                  # barcode -> cell_type
    counts_file_path='counts_normalized.h5ad',  # normalized, NOT scaled
    counts_data='hgnc_symbol',
    threshold=0.1, iterations=1000, pvalue=0.05,
    score_interactions=True, threads=4, output_path='cpdb_out')
# DEG-driven escape from one-vs-rest: cpdb_degs_analysis_method.call(..., degs_file_path=...)
```

## Pathway Probability (CellChat v2)

**Goal:** Summarize communication at the signaling-pathway level with sender/receiver roles.

**Approach:** Build the object, pick a database subset, identify over-expressed interactions, compute the mass-action probability with trimean, filter tiny populations, aggregate to pathways, then compute centrality for role analysis. Order matters.

```r
library(CellChat)

cellchat <- createCellChat(object = seurat_obj, group.by = 'cell_type')
cellchat@DB <- CellChatDB.human   # or CellChatDB.mouse; subsetDB(..., search='Secreted Signaling') to restrict
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = 'triMean')   # trimean ~25% truncated mean: conservative
cellchat <- filterCommunication(cellchat, min.cells = 10)   # drop populations under 10 cells
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = 'netP')   # sender/receiver/mediator roles
# Viz: netVisual_aggregate(signaling='WNT'), netVisual_bubble(), netAnalysis_signalingRole_heatmap()
```

## Downstream Ligand-Activity (NicheNet)

**Goal:** Identify which sender ligand best explains the receiver's observed transcriptional response - the distinct, better-grounded question.

**Approach:** Define a receiver gene set of interest (DE genes from a condition contrast), restrict to ligands expressed in senders with receptors expressed in the receiver, and rank ligands by how well their predicted regulatory targets recover that gene set (AUPR).

```r
library(nichenetr)
library(Seurat)
library(tidyverse)

ligand_target_matrix <- readRDS('ligand_target_matrix.rds')
lr_network <- readRDS('lr_network.rds')

# Receiver gene set: garbage in -> garbage out; a noisy/batch-confounded DE list invalidates the ranking
geneset_oi <- FindMarkers(seurat_obj, ident.1 = 'activated_T', ident.2 = 'naive_T') %>%
    filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>% rownames()
background <- get_expressed_genes('T_cell', seurat_obj, pct = 0.10)

expressed_ligands <- intersect(unique(lr_network$from), get_expressed_genes(c('Macrophage', 'Dendritic'), seurat_obj, 0.10))
expressed_receptors <- intersect(unique(lr_network$to), background)
potential_ligands <- lr_network %>% filter(from %in% expressed_ligands, to %in% expressed_receptors) %>% pull(from) %>% unique()

ligand_activities <- predict_ligand_activities(
    geneset = geneset_oi, background_expressed_genes = background,
    ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

# Current model ranks by aupr_corrected (AUPR is the headline metric; v1 used pearson)
best_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)
```

## Threshold and Permutation Rationale

| Parameter | Default | Rationale |
|-----------|---------|-----------|
| `expr_prop` / `threshold` | 0.10 | A gene expressed in <10% of a cluster is mostly dropout; below this, scores are noise - but real low-abundance signaling is also discarded (the "not necessary" side of the proxy) |
| `n_perms` / `iterations` | 1000 | Stable label-permutation p-values; 100 is fine for exploration, 1000 for reporting; the p-value is about label shuffling, not binding |
| `min.cells` (CellChat) | 10 | Populations under ~10 cells give unstable mean expression and inflated probabilities |
| trimean (CellChat) | type='triMean' | 25% truncated mean is conservative, yielding fewer, higher-confidence calls than CellPhoneDB's mean |
| `aupr_corrected` top-N | 30 | NicheNet ligand cutoff is a display choice, not a significance threshold; inspect the activity-score elbow |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `KeyError: 'liana_rank'` | `rank_aggregate` outputs `magnitude_rank` and `specificity_rank`, not a single combined rank | Filter on `specificity_rank` and/or `magnitude_rank` |
| One dominant cluster is the hub of every network | Abundance drives the permutation null; ambient RNA inflates its secreted ligands | Decontaminate ambient RNA, down-sample, check counts vs cluster size |
| Findings flip when the database changes | Resource choice moves results as much as method (Dimitrov 2022) | Report the resource and show the key pair survives >=2 resources |
| Contact-dependent pair (Notch-DLL, ephrin) called between non-adjacent types | Membrane-bound ligands scored as if secreted; no geometry in dissociated data | Restrict to secreted signaling or use a spatial method with proximity |
| AP-1 / heat-shock "stress signaling" everywhere | Dissociation-induced FOS/JUN/HSPA modules treated as ligands | Flag/regress the stress module before scoring |
| NicheNet ligand ranking looks random | Receiver gene set is noisy or batch-confounded; or pathway is inactive in that lineage (static prior) | Clean the DE contrast; treat top ligands as "consistent with the response under a generic prior" |
| Mouse CellPhoneDB run returns almost nothing | CellPhoneDB DB is human-only | Map orthologs or use CellChatDB.mouse / LIANA `mouseconsensus` |
| More interactions claimed in condition B than A | Interaction counts scale with cell number and depth | Compare score magnitudes or use CellChat differential / Tensor-cell2cell, not raw counts |

## Related Skills

- single-cell/cell-annotation - Cell-type labels define senders and receivers; annotation resolution is a hidden CCC hyperparameter
- single-cell/clustering - Cluster granularity changes who is "specific"; fix it before running CCC
- single-cell/doublet-detection - Doublets create fake co-expressing cells that masquerade as senders-receivers
- single-cell/preprocessing - Ambient-RNA decontamination and stress-gene handling happen here, before CCC
- single-cell/metabolite-communication - Metabolite-mediated CCC (enzyme-sensor) as the doubly-inferred counterpart to ligand-receptor
- spatial-transcriptomics/spatial-communication - Proximity-constrained CCC when spatial coordinates are available
- pathway-analysis/go-enrichment - Functional enrichment of NicheNet target genes or interacting receptors
- differential-expression/deseq2-basics - Pseudobulk DE to build the receiver gene set NicheNet requires

## References

- Vento-Tormo R, Efremova M, et al. Single-cell reconstruction of the early maternal-fetal interface in humans. Nature 563:347-353 (2018). [original CellPhoneDB]
- Efremova M, Vento-Tormo M, Teichmann SA, Vento-Tormo R. CellPhoneDB: inferring cell-cell communication from combined expression of multi-subunit ligand-receptor complexes. Nat Protoc 15:1484-1506 (2020). [statistical method]
- Jin S, et al. Inference and analysis of cell-cell communication using CellChat. Nat Commun 12:1088 (2021).
- Browaeys R, Saelens W, Saeys Y. NicheNet: modeling intercellular communication by linking ligands to target genes. Nat Methods 17(2):159-162 (2020).
- Dimitrov D, et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nat Commun 13:3224 (2022). [LIANA, discordance]
- Dimitrov D, et al. LIANA+ provides an all-in-one framework for cell-cell communication inference. Nat Cell Biol 26:1613-1622 (2024).
- Hou R, et al. Predicting cell-to-cell communication networks using NATMI. Nat Commun 11:5011 (2020).
- Cang Z, Nie Q, et al. Screening cell-cell communication in spatial transcriptomics via collective optimal transport [COMMOT]. Nat Methods 20:218-228 (2023).
- Palla G, et al. Squidpy: a scalable framework for spatial omics analysis. Nat Methods 19:171-178 (2022).
- Luo J, et al. ESICCC: evaluation, selection, and integration of cell-cell communication inference methods. Genome Res 33(10):1788-1805 (2023). [benchmark]
- Young MD, Behjati S. SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data. GigaScience 9(12):giaa151 (2020).
- van den Brink SC, et al. Single-cell sequencing reveals dissociation-induced gene expression in tissue subpopulations. Nat Methods 14(10):935-936 (2017).
