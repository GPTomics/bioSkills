---
name: bio-single-cell-perturb-seq
description: Analyze Perturb-seq / CROP-seq single-cell CRISPR screens. Use when assigning guides as a mixture problem, removing non-perturbed escaper cells with Mixscape, choosing a calibrated test (SCEPTRE conditional resampling) over naive DE, quantifying effect size with E-distance, separating compositional shifts from within-state expression change, or judging whether a perturbation-prediction foundation model actually beats a baseline.
tool_type: python
primary_tool: Pertpy
---

## Version Compatibility

Reference examples tested with: pertpy 0.9+, scanpy 1.10+, anndata 0.10+, sceptre 0.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Perturb-seq Analysis

**"Analyze my Perturb-seq CRISPR screen"** -> Assign guides, remove cells that received a guide but were not perturbed, test each perturbation with a calibrated method, and separate "moves cells" from "changes cells".
- Python: `pertpy.pp.GuideAssignment`, `pertpy.tl.Mixscape`, `pertpy.tl.Distance`/`DistanceTest`, `pertpy.tl.Milo`/`Sccoda`
- R: `sceptre` (conditional-resampling test), `Seurat` Mixscape, `scMAGeCK`

## Governing Principle

Guide assignment is a mixture problem, not a threshold. Each cell's per-guide UMI vector mixes true integration with ambient guide contamination (free transcripts, index hopping, doublets), and the ambient pool is structured: it is dominated by whichever guides are most abundant in the library, so a flat UMI cutoff preferentially mis-assigns cells to common guides and calls rare-guide cells negative. Call guides by a per-guide background/foreground mixture posterior, and report the perturbed fraction. MOI changes the meaning: low-MOI (~1 guide/cell) gives clean single-gene attribution but discards 70-90% of cells; high-MOI is for combinatorial designs but measures every single-gene effect in a co-perturbed background.

Assignment is not effective perturbation. A cell can carry a guide yet be transcriptionally wild-type: incomplete CRISPR-KO editing, in-frame indels, escapers, or weak CRISPRi knockdown. The "perturbed" population is a mixture of truly perturbed and effectively-wild-type cells, which attenuates every effect-size estimate toward the null. Mixscape removes the non-perturbed cells via a local non-targeting-neighbor perturbation signature before testing. Deep caveat: an all-NP result is not evidence the gene is non-functional, because it is confounded with low guide efficiency; Mixscape cannot distinguish "no phenotype" from "no editing".

Naive DE is miscalibrated by depth and pseudoreplication. The probability of detecting a guide covaries with sequencing depth, and depth also drives expression detection, so a plain Wilcoxon/NB test between guide-positive and NT cells has inflated type-I error. SCEPTRE fixes this with conditional resampling: it models P(cell receives this guide | technical covariates incl. depth) and resamples the assignment to build a calibrated null, robust to misspecification of the expression model. Separately, treating thousands of cells from one transfection as independent replicates inflates significance (pseudoreplication, Squair 2021); the replication unit is the transfection, so pseudobulk-per-replicate is required for calibrated inference.

E-distance is the modern effect-size. Energy distance in PCA space, `E = 2*sigma_between - sigma_within_X - sigma_within_Y`, measures separation magnitude (not direction or mechanism), with a permutation E-test. It is interpretable only relative to a fixed embedding and is not comparable across studies with different pipelines.

Separate the compositional shift from the within-state change. A perturbation can (a) shift the proportions of pre-existing cell states (differential abundance) without changing any state's program, or (b) change expression within a state (differential expression). A perturbation that only redistributes cells produces a large pseudobulk "DE signature" that is entirely a composition artifact. These need different tools and answer different questions; report both.

## Guide Assignment: Mixture vs Threshold

| Method | Model | Use when | Fails when |
|---|---|---|---|
| Mixture (posterior) | Per-guide 2-component mixture (background Poisson + foreground Gaussian); `pt.pp.GuideAssignment.assign_mixture_model` | Default; ambient varies by guide abundance; low and high MOI | Very few cells per guide (mixture unstable); verify against NT contamination floor |
| Threshold | Flat UMI cutoff; `assign_by_threshold` | Quick sanity check, uniform high-signal libraries | Ambient scales with abundant guides -> mis-assigns to common guides, calls rare-guide cells negative |

Cell Ranger and Replogle's `guide_calling` fit mixtures on log counts; require a minimum dominant-guide UMI fraction, not just an absolute count, and gate doublets (they masquerade as combinatorial cells).

## Method Decision Table (Testing and Effect Size)

| Method | What it answers | Use when | Fails when |
|---|---|---|---|
| Mixscape (pertpy/Seurat) | Which cells were effectively perturbed; per-perturbation DE after removing escapers | CRISPR-KO with heterogeneous editing; need escaper removal | All-NP confounded with low guide efficiency; KO posteriors not comparable across targets |
| SCEPTRE | Calibrated perturbation-gene association | Rigorous testing under the depth confounder; element-level screens | Needs the assignment model roughly right; conservative by design |
| scMAGeCK (LR / RRA) | Per-gene effect across many genes; high-MOI deconvolution | Multi-guide cells; ridge-regression effect estimates | NEGCTRL choice defines the null; runs on scale.data so covariates propagate |
| E-distance + E-test (pertpy) | Effect-size magnitude; perturbation similarity | Ranking/clustering perturbations by how far they move cells | Embedding-dependent, not cross-study comparable; floored by permutation count |
| Pseudobulk DE (DESeq2/edgeR) | Average within-state program change | >=2-3 biological replicates per condition | One replicate per guide -> no valid inference; sum raw counts, not means |
| Milo / scCODA / Augur | Differential abundance / composition | "Does the perturbation move cells across states?" | Conflated with within-state DE if reported alone |

Verify the current best-practice default and parameter names against the installed pertpy/sceptre docs before committing; the APIs drift across releases.

## Foundation-Model Reality Check

This is settled as of 2026, not hype. scGPT, Geneformer, scFoundation, scBERT, UCE are pretrained with masked-expression objectives and learn the co-expression manifold of observational data; perturbation prediction is a causal/interventional question, and there is no theorem that co-expression transfers to intervention. The empirical result across benchmarks: none reliably beat trivial baselines on unseen perturbations (Ahlmann-Eltze 2025; Kernfeld 2025; Csendes 2025).
- For unseen single perturbations, predicting the mean perturbed profile across training perturbations is hard to beat; for unseen doubles, an additive model (sum the two single effects) captures most variance because genetic interactions are the exception.
- All-gene MSE/correlation is dominated by unchanged genes, so a "predict no change"/mean model scores deceptively high. Evaluate on DE genes against mean/additive baselines.
- Train-test leakage inflates reported success: random cell-level splits put the same perturbation in train and test. True generalization holds out entire perturbations (and ideally cell-type contexts), not random cells.

The defensible reviewer stance: demand whole-perturbation holdout, DE-gene metrics, and an explicit additive/mean baseline. Without these, a positive result is not credible.

## Guide Assignment (pertpy)

**Goal:** Call which guide each cell actually received using a mixture posterior, not a flat threshold.

**Approach:** Fit a per-guide Poisson-Gaussian mixture to the guide-count modality and assign by posterior, allowing negative and multi-guide calls.

```python
import pertpy as pt
import scanpy as sc

gdo = mdata.mod['gdo']                       # guide-count modality (cells x guides)
gdo.layers['counts'] = gdo.X.copy()

ga = pt.pp.GuideAssignment()
ga.assign_mixture_model(gdo, assigned_guides_key='assigned_guide')   # background Poisson + foreground Gaussian
# Inspect NT/abundant-guide UMI distributions as a contamination floor before trusting calls
ga.plot_heatmap(gdo, layer='counts')
```

## Mixscape: Remove Non-Perturbed Escapers (pertpy)

**Goal:** Separate effectively-perturbed (KO) from non-perturbed (NP) cells before any DE.

**Approach:** Build a local perturbation signature by subtracting each cell's NT neighbors, then fit a per-target 2-component mixture to classify cells; drop NP cells.

```python
ms = pt.tl.Mixscape()
ms.perturbation_signature(adata, pert_key='perturbation', control='NT', n_neighbors=20)   # pert_key here = the broad perturbed-vs-control column
ms.mixscape(adata, pert_key='target_gene', control='NT', layer='X_pert')   # pert_key here = the per-target column (intentionally different); renamed from labels; writes adata.obs['mixscape_class_global'] KO/NP/NT
# An all-NP target is confounded with low guide efficiency: report perturbed fraction, do not call the gene non-functional
adata.obs['mixscape_class_global'].value_counts()
```

## E-distance and the E-test (pertpy)

**Goal:** Quantify how far each perturbation moves cells and test it against a permutation null.

**Approach:** Compute energy distance in a fixed PCA embedding; pin the embedding and metric, and run the permutation E-test against the control.

```python
sc.pp.pca(adata, n_comps=50)
dist = pt.tl.Distance(metric='edistance', obsm_key='X_pca')   # pin obsm; sqeuclidean vs euclidean default changed across versions
pairwise = dist.pairwise(adata, groupby='target_gene')

etest = pt.tl.DistanceTest('edistance', n_perms=1000)         # smallest p ~ 1/(n_perms+1); crushed by multiple testing
results = etest(adata, groupby='target_gene', contrast='NT')
```

## SCEPTRE: Calibrated Testing (R)

**Goal:** Test perturbation-gene associations with calibration verified on the data itself.

**Approach:** Import counts and guide matrices, set parameters, assign guides by mixture, then run the calibration check (negative controls) before the discovery analysis.

```r
library(sceptre)

obj <- import_data(response_matrix = rna_counts, grna_matrix = grna_counts,
                   grna_target_data_frame = grna_targets, moi = 'low')
obj <- set_analysis_parameters(obj, discovery_pairs = pairs)
obj <- assign_grnas(obj, method = 'mixture')        # mixture | thresholding | maximum
obj <- run_qc(obj)
obj <- run_calibration_check(obj)                   # negative-control pairs must be well-calibrated FIRST
obj <- run_discovery_analysis(obj)
results <- get_result(obj, analysis = 'discovery_analysis')
```

## Pseudobulk DE (Within-State Change)

**Goal:** Test the average program change per perturbation with valid biological replication.

**Approach:** Sum RAW counts per (target gene, replicate), filter tiny pseudobulk samples, then run DESeq2/edgeR; this respects the replication unit and avoids pseudoreplication.

```python
import pertpy as pt

adata.layers['counts'] = adata.layers.get('counts', adata.X.copy())   # stash RAW counts before any log1p
pb = pt.tl.PseudobulkSpace()
pdata = pb.compute(adata, target_col='target_gene', groups_col='replicate', layer_key='counts', mode='sum')   # sum RAW counts, not .X (log-normalized)
# Drop pseudobulk samples below ~10 cells (verify the per-sample cell-count obs column name with help(pb.compute))
# Hand pdata to pertpy EdgeR / pydeseq2 with design ~ replicate + target_gene; needs >=2-3 replicates per condition
```

One transfection per guide means no valid biological-replicate inference exists; using guides targeting the same gene as pseudo-replicates partially helps but conflates guide-specific off-targets.

## Compositional vs Expression (Separate the Questions)

**Goal:** Decide whether a perturbation moves cells across states or changes a state's program.

**Approach:** Run a differential-abundance test (Milo neighborhoods or scCODA) for composition, and report it alongside the within-state pseudobulk DE.

```python
milo = pt.tl.Milo()
mdata_milo = milo.load(adata)
milo.make_nhoods(mdata_milo['rna'])
milo.count_nhoods(mdata_milo, sample_col='replicate')
milo.da_nhoods(mdata_milo, design='~ target_gene')   # differential abundance: does the perturbation shift proportions?
```

## Common Errors

| Symptom | Cause | Fix |
|---|---|---|
| Rare-guide cells called negative | Flat UMI threshold; ambient biased to abundant guides | Mixture-model assignment by posterior; require a dominant-guide UMI fraction |
| Effect sizes weaker than expected | Escapers/incomplete KO dilute the perturbed population | Run Mixscape, remove NP cells, report perturbed fraction |
| "Gene is non-functional" from all-NP | All-NP confounds no-phenotype with no-editing | Do not claim non-functional; check guide efficiency independently |
| Hundreds of "significant" hits | Naive Wilcoxon/NB miscalibrated by depth + pseudoreplication | SCEPTRE conditional resampling; pseudobulk-per-replicate DE |
| Huge DE signature but no program change | Perturbation only redistributes cells across states | Run Milo/scCODA; attribute the signal to composition |
| E-distances disagree with another paper | Embedding/metric/PC count differ; default metric changed | Pin pertpy version, obsm key, and cell_wise_metric; do not cross-compare |
| Combinatorial cells everywhere | Doublets masquerade as multi-guide | Gate doublets (Scrublet/scDblFinder) before multi-guide analysis |
| Foundation model "beats" baselines | Cell-level split leakage; all-gene metric hides failure | Hold out whole perturbations; score DE genes vs additive/mean baseline |

## Related Skills

- single-cell/preprocessing - scRNA-seq QC and normalization upstream of the screen
- single-cell/doublet-detection - gating doublets before multi-guide analysis
- single-cell/markers-annotation - interpreting per-perturbation DE genes
- single-cell/differential-abundance - compositional shift testing (Milo/scCODA) for perturbations that change cell-state proportions
- single-cell/batch-integration - multi-sample/replicate integration
- crispr-screens/mageck-analysis - bulk CRISPR screen analysis (MAGeCK RRA/MLE)
- crispr-screens/perturb-seq-analysis - related single-cell CRISPR screen workflow
- differential-expression/deseq2-basics - pseudobulk DESeq2 testing on summed counts
- pathway-analysis/go-enrichment - pathway interpretation of perturbation signatures

## References

Dixit A, Parnas O, Li B, et al. Perturb-Seq: dissecting molecular circuits with scalable single-cell RNA profiling of pooled genetic screens. Cell 167(7):1853-1866 (2016).
Datlinger P, Rendeiro AF, Schmidl C, et al. Pooled CRISPR screening with single-cell transcriptome readout (CROP-seq). Nat Methods 14(3):297-301 (2017).
Replogle JM, Norman TM, Xu A, et al. Combinatorial single-cell CRISPR screens by direct guide RNA capture and targeted sequencing. Nat Biotechnol 38(8):954-961 (2020).
Papalexi E, Mimitou EP, Butler AW, et al. Characterizing the molecular regulation of inhibitory immune checkpoints with multimodal single-cell screens (Mixscape). Nat Genet 53(3):322-331 (2021).
Yang L, Zhu Y, Yu H, et al. scMAGeCK links genotypes with multiple phenotypes in single-cell CRISPR screens. Genome Biol 21:19 (2020).
Barry T, Wang X, Morris JA, Roeder K, Katsevich E. SCEPTRE improves calibration and sensitivity in single-cell CRISPR screen analysis. Genome Biol 22:344 (2021).
Squair JW, Gautier M, Kathe C, et al. Confronting false discoveries in single-cell differential expression. Nat Commun 12:5692 (2021).
Peidli S, Green TD, Shen C, et al. scPerturb: harmonized single-cell perturbation data (E-distance). Nat Methods 21(3):531-540 (2024).
Heumos L, Ji Y, May L, et al. Pertpy: an end-to-end framework for perturbation analysis. Nat Methods 23(2):350-359 (2026).
Dann E, Henderson NC, Teichmann SA, Morgan MD, Marioni JC. Differential abundance testing on single-cell data using k-nearest neighbor graphs (Milo). Nat Biotechnol 40(2):245-253 (2022).
Ahlmann-Eltze C, Huber W, Anders S. Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines. Nat Methods 22(8):1657-1661 (2025).
Kernfeld E, Yang Y, Weinstock JS, et al. A comparison of computational methods for expression forecasting. Genome Biol 26:388 (2025).
Csendes G, et al. Benchmarking foundation cell models for post-perturbation RNA-seq prediction. BMC Genomics 26:393 (2025).
