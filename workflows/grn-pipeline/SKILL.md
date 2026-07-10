---
name: bio-workflows-grn-pipeline
description: Orchestrates gene regulatory network inference from processed single-cell data to regulons and in-silico perturbation, via pySCENIC (RNA-only GRNBoost2 -> cisTarget -> AUCell), SCENIC+ (multiome cisTopic -> pycistarget -> eGRN), and CellOracle perturbation. Use when recognizing that an inferred GRN is UNDIRECTED by default and reporting only the evidence tier delivered (co-expression vs motif-pruned vs enhancer-resolved vs perturbation), matching species/assembly/namespace across the TF-list + cisTarget DB + motif2TF annotation, feeding RAW counts of the cleaned/doublet-free/batch-controlled cells (never imputed/batch-corrected values), running the cisTarget pruning that buys directionality (modules are not regulons without it), or choosing the RNA-only vs multiome path. Hands mechanism to the gene-regulatory-networks component skills; not a re-teach of any single step.
tool_type: python
primary_tool: pySCENIC
workflow: true
depends_on:
  - gene-regulatory-networks/scenic-regulons
  - gene-regulatory-networks/multiomics-grn
  - gene-regulatory-networks/perturbation-simulation
  - single-cell/clustering
qc_checkpoints:
  - after_grn_inference: "50-500 regulons detected, known TFs present"
  - after_activity_scoring: "AUCell scores separate known cell types"
  - after_perturbation: "Predicted shifts match known biology"
---

## Version Compatibility

Reference examples tested with: pySCENIC 0.12+, arboreto 0.1.6+, ctxcore 0.2+, pycisTopic 2.0+, pycistarget 1.0+, SCENIC+ 1.0a1 (Snakemake CLI), CellOracle 0.18+, anndata 0.10+, pandas 2.2+, scanpy 1.10+, scipy 1.12+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: SCENIC+ is now a Snakemake pipeline (`scenicplus init_snakemake`); the pre-2024 manual `create_SCENICPLUS_object`/`build_grn` object API is deprecated. GRNBoost2's arboreto dask backend is the #1 operational landmine — use the bundled multiprocessing if the dask cluster hangs. The TF-list, cisTarget ranking DB, and motif2TF `.tbl` must all be the SAME species + assembly + collection vintage. Confirm in-tool before quoting.

# Gene Regulatory Network Pipeline

**"Infer gene regulatory networks from my single-cell data"** -> Orchestrate pySCENIC regulon inference (GRNBoost2, cisTarget, AUCell), CellOracle perturbation simulation, and regulon-based cell type characterization.

This is a workflow skill: it owns the chaining decisions and hand-offs, not the internals of any one step.

## The governing principle

1. **An inferred GRN is an UNDIRECTED association graph by default; directionality is IMPORTED, and only the evidence tier actually delivered should be reported.** GRNBoost2 co-expression is tier-1 (undirected); the cisTarget motif-pruning step (Path A step 2) is what buys directionality and discards indirect edges — modules before ctx are NOT regulons, calling them so is a category error. SCENIC+ adds enhancer resolution; perturbation adds causal direction. Do not write tier-6 "master regulator drives X" prose over a tier-1 co-expression result.
2. **Species/assembly/namespace must match across the TF-list, the cisTarget ranking DB, and the motif2TF `.tbl` — all three.** A mismatch (mouse genes in an hg38 DB; feather v1 DB with a v10 motif annotation) yields near-empty regulons. Also commit and report the search-space window (500bp/100bp proximal vs TSS±10kb) — results are not comparable across windows.
3. **Feed RAW counts of the CLEANED, doublet-free, batch-controlled cells — never imputed or batch-corrected values.** GRNBoost2 on imputed counts inflates correlations (imputation smooths neighbors into agreement); on batch-corrected values a batch module can pass motif enrichment by chance; doublets create a fake "hybrid regulator." Run SCENIC ONCE on the integrated object.
4. **Validate a regulon with an ORTHOGONAL modality, not the TF's own mRNA.** AUCell activity != TF expression; correlating a regulon's activity with its TF's expression is circular (and activity is dropout-robust while the TF mRNA may read zero). GRNBoost2 is stochastic — run multiple seeds and keep recurrent links (Van de Sande 2020).

## Made-once commitments

| Commitment | Consequence inherited downstream |
|------------|----------------------------------|
| Species + assembly + gene namespace (HGNC vs MGI; hg38 vs mm10) | TF list, cisTarget DB, motif2TF `.tbl` must all match, or regulons are near-empty |
| cisTarget DB vintage + search window (proximal vs TSS±10kb; gene- vs region-based) | Which edges survive ctx pruning; results not comparable across windows/vintages |
| Input matrix identity (RAW counts, from the cleaned/doublet-free/integrated object) | Every adjacency and regulon; imputed/batch values fabricate edges |
| RNA-only vs multiome availability | Which path is possible: SCENIC+ REQUIRES paired multiome; RNA-only -> pySCENIC or CellOracle-with-prebuilt-base-GRN |

## Pipeline Overview

```
Processed AnnData (QC'd, normalized, clustered)
    |
    +----- RNA only? -------> Path A: pySCENIC (3-step)
    |                              |
    |                              v
    |                         [1. GRNBoost2] ----> TF-target adjacencies
    |                              |
    |                              v
    |                         [2. RcisTarget] ---> Regulon pruning (motif enrichment)
    |                              |
    |                              v
    |                         [3. AUCell] -------> Regulon activity scoring
    |
    +----- Multiome? -------> Path B: SCENIC+
    |                              |
    |                              v
    |                         [1. cisTopic] -----> Topic modeling on ATAC
    |                              |
    |                              v
    |                         [2. pycistarget] --> Enhancer-TF mapping
    |                              |
    |                              v
    |                         [3. SCENIC+] ------> eGRN construction
    |
    +---> [CellOracle Perturbation Simulation] (either path)
              |
              v
         Perturbation scores + predicted cell state shifts
```

## Path A: pySCENIC (RNA-Only)

### Step 1: GRN Inference with GRNBoost2

```python
import scanpy as sc
import pandas as pd
from arboreto.algo import grnboost2

adata = sc.read_h5ad('processed.h5ad')

# Extract expression matrix (raw counts recommended for GRNBoost2)
expr_matrix = pd.DataFrame(
    adata.raw.X.toarray() if hasattr(adata.raw.X, 'toarray') else adata.raw.X,
    index=adata.obs_names, columns=adata.raw.var_names
)

# TF list from cisTarget resources
# Human: https://resources.aertslab.org/cistarget/tf_lists/
tf_names = pd.read_csv('allTFs_hg38.txt', header=None)[0].tolist()
tf_names = [tf for tf in tf_names if tf in expr_matrix.columns]

adjacencies = grnboost2(expr_matrix, tf_names=tf_names, seed=42, verbose=True)
adjacencies.to_csv('adjacencies.tsv', sep='\t', index=False, header=False)
```

### Step 2: Regulon Pruning with RcisTarget

```python
from pyscenic.prune import prune2df, df2regulons
from pyscenic.utils import modules_from_adjacencies
from ctxcore.rnkdb import FeatherRankingDatabase

# cisTarget databases (~10 GB each, download once)
# Human: hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
# Mouse: mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
# FeatherRankingDatabase(fname, name) -- `name` is REQUIRED (no default) in ctxcore
dbs = [FeatherRankingDatabase(db, name=os.path.splitext(os.path.basename(db))[0]) for db in [
    'hg38_500bp_up_100bp_down.genes_vs_motifs.rankings.feather',
    'hg38_10kbp_up_10kbp_down.genes_vs_motifs.rankings.feather'
]]

motif_annotations = 'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl'

# Build co-expression modules as Regulon objects. prune2df reads module.transcription_factor,
# which a bare GeneSignature lacks (AttributeError); modules_from_adjacencies applies pySCENIC's
# standard top-target/importance thresholds and returns the Regulon objects prune2df expects.
modules = list(modules_from_adjacencies(adjacencies, expr_matrix))

# Prune modules using motif enrichment
# NES threshold 3.0 (default); rank_threshold=5000 matches the CLI (prune2df default is 1500).
df_motifs = prune2df(dbs, modules, motif_annotations, rank_threshold=5000, num_workers=8)
regulons = df2regulons(df_motifs)

print(f'Discovered {len(regulons)} regulons')
```

### Step 3: AUCell Activity Scoring

```python
from pyscenic.aucell import aucell

auc_matrix = aucell(expr_matrix, regulons, num_workers=8)

adata.obsm['X_aucell'] = auc_matrix.loc[adata.obs_names].values
adata.uns['regulon_names'] = [r.name for r in regulons]
```

### QC Checkpoint: GRN Inference

```python
def validate_grn(regulons, auc_matrix, adata, cell_type_key='cell_type'):
    '''
    QC gates after GRN inference.
    - 50-500 regulons is typical range
    - Known lineage TFs should appear (e.g., PAX6 in neurons, GATA1 in erythroid)
    - AUCell scores should separate known cell types
    '''
    n_regulons = len(regulons)
    regulon_names = [r.name for r in regulons]

    # Gate 1: Regulon count
    if n_regulons < 50:
        print(f'WARNING: Only {n_regulons} regulons. Check TF list or lower NES threshold.')
    elif n_regulons > 500:
        print(f'WARNING: {n_regulons} regulons found. Consider stricter pruning.')
    else:
        print(f'OK: {n_regulons} regulons in expected range (50-500)')

    # Gate 2: Known TFs present
    known_tfs = ['PAX6', 'SOX2', 'GATA1', 'SPI1', 'FOXP3', 'TBX21', 'EBF1']
    # df2regulons names regulons 'PAX6(+)' / 'PAX6(-)'; strip the suffix or this gate never fires
    regulon_bases = {name.split('(')[0] for name in regulon_names}
    found = [tf for tf in known_tfs if tf in regulon_bases]
    print(f'Known lineage TFs found: {found}')

    # Gate 3: AUCell separates cell types
    import scipy.stats as stats
    cell_types = adata.obs[cell_type_key].unique()
    if len(cell_types) >= 2:
        ct1_idx = adata.obs[cell_type_key] == cell_types[0]
        ct2_idx = adata.obs[cell_type_key] == cell_types[1]
        n_differential = 0
        for i, rname in enumerate(regulon_names[:min(50, len(regulon_names))]):
            stat, pval = stats.mannwhitneyu(
                auc_matrix.values[ct1_idx, i], auc_matrix.values[ct2_idx, i]
            )
            if pval < 0.01:
                n_differential += 1
        print(f'Differentially active regulons between top 2 types: {n_differential}/50')

    return n_regulons
```

## Path B: SCENIC+ (Multiome)

### Step 1: ATAC Topic Modeling with cisTopic

```python
import pycisTopic
from pycisTopic.cistopic_class import create_cistopic_object
from pycisTopic.lda_models import run_cgs_models

# Create cisTopic object from fragments
cistopic_obj = create_cistopic_object(
    fragment_matrix=adata_atac.X.T,   # cisTopic wants regions x cells; AnnData .X is cells x regions
    cell_names=adata_atac.obs_names.tolist(),
    region_names=adata_atac.var_names.tolist()
)

# Run LDA topic modeling
# n_topics: test range around expected cell types (e.g., 2x number of clusters)
models = run_cgs_models(
    cistopic_obj,
    n_topics=[10, 20, 30, 40, 50],
    n_cpu=8, n_iter=300, random_state=42
)

# evaluate_models plots the model-selection metrics; read the elbow and pass that topic COUNT
# as select_model (an int, not True -- True==1 would select a non-existent 1-topic model).
from pycisTopic.lda_models import evaluate_models
model = evaluate_models(models, select_model=40, return_model=True)
cistopic_obj.add_LDA_model(model)
```

### Step 2: Enhancer-TF Mapping

```python
import pyranges as pr
from pycistarget.utils import region_names_to_coordinates
from pycistarget.motif_enrichment_cistarget import run_cistarget
from pycisTopic.topic_binarization import binarize_topics

region_bin = binarize_topics(cistopic_obj, method='otsu')   # dict of DataFrames keyed by topic (region names in the index)

# run_cistarget needs a dict of pyranges.PyRanges, not the raw binarized DataFrames.
region_sets = {topic: pr.PyRanges(region_names_to_coordinates(region_bin[topic].index.tolist()))
               for topic in region_bin}

# Run motif enrichment on accessible regions. The first arg is the cisTarget ranking DB: pass the
# feather path and run_cistarget instantiates cisTargetDatabase itself. Prebuilt DBs at
# https://resources.aertslab.org/cistarget/ . The parameter is spelled `specie`, not `species`.
CTX_DB = '/path/to/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather'
cistarget_results = run_cistarget(
    CTX_DB,
    region_sets=region_sets,
    specie='homo_sapiens',
    auc_threshold=0.005,
    nes_threshold=3.0,
    rank_threshold=0.05,
    n_cpu=8
)
```

### Step 3: eGRN Construction

**Goal:** Assemble eRegulons (TF -> enhancer -> gene triplets) from the multiome data.

**Approach:** Current SCENIC+ runs topic modeling, motif enrichment, and eGRN construction through one Snakemake pipeline; the deprecated manual `create_SCENICPLUS_object`/`build_grn` API (and pre-2024 tutorials) should not be used. See gene-regulatory-networks/multiomics-grn for the full pipeline and the peak-to-gene caveats.

```bash
# Scaffold, edit the config (point at fragments, scRNA AnnData, cell-type labels, databases),
# then run from inside the Snakemake directory.
scenicplus init_snakemake --out_dir scenicplus_run
# edit scenicplus_run/Snakemake/config/config.yaml
cd scenicplus_run/Snakemake && snakemake --cores 16
```

```python
# Read the resulting direct (high-confidence) eRegulon table (filename is config-/version-
# dependent, so resolve it by glob).
import glob, pandas as pd
eregulons = pd.read_csv(glob.glob('scenicplus_run/**/eRegulon*direct*.tsv', recursive=True)[0], sep='\t')
print(f'eRegulons: {eregulons["TF"].nunique()} enhancer-driven regulators')
```

## CellOracle Perturbation Simulation

**Goal:** Predict the direction cells move under a TF knockout, as a hypothesis (direction, not calibrated magnitude).

**Approach:** CellOracle needs a base GRN (a TF-target scaffold from motif scanning of accessible regions, not the pySCENIC adjacencies), then learns per-cluster weights, propagates a forced expression shift, and projects it onto the cell-state graph. See gene-regulatory-networks/perturbation-simulation for the base-GRN construction and the local-linear / direction-only caveats.

```python
import celloracle as co
import numpy as np

oracle = co.Oracle()
oracle.import_anndata_as_raw_count(adata=adata, cluster_column_name='cell_type',
                                   embedding_name='X_umap')

# Base GRN = motif-scanned accessible regions (preferred) or a prebuilt CellOracle base GRN;
# this is NOT the pySCENIC adjacencies. See multiomics-grn / perturbation-simulation.
base_grn = co.data.load_human_promoter_base_GRN()   # `version` must match the genome build
oracle.import_TF_data(TF_info_matrix=base_grn)

oracle.perform_PCA()
k = int(0.025 * oracle.adata.n_obs)
oracle.knn_imputation(n_pca_dims=50, k=k, balanced=True, b_sight=k * 8, b_maxl=k * 4)

# Learn context-specific weights, then fit the simulation GRN.
links = oracle.get_links(cluster_name_for_GRN_unit='cell_type', alpha=10)
links.filter_links(p=0.001, weight='coef_abs', threshold_number=2000)
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)

# Simulate TF knockout (0.0) and project the shift onto the embedding.
oracle.simulate_shift(perturb_condition={'MYC': 0.0}, n_propagation=3)
oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
oracle.calculate_embedding_shift(sigma_corr=0.05)
shift = np.sqrt((oracle.delta_embedding ** 2).sum(axis=1))
```

### QC Checkpoint: Perturbation

```python
def validate_perturbation(oracle, perturbed_tf, expected_affected_cluster=None):
    '''
    QC gate: perturbation shifts should match known biology.
    - Transition probabilities should show directional shift
    - If expected_affected_cluster known, check it shows largest change
    '''
    import numpy as np, pandas as pd
    # Shift magnitude per cell from the simulated embedding shift (delta_embedding).
    shift = np.sqrt((oracle.delta_embedding ** 2).sum(axis=1))
    # observed=True: cell_type is categorical, and the default retains filtered-out categories as NaN rows.
    # Sort here, not at print time: the gate below reads index[:3], which is category order until sorted.
    mean_shift = pd.Series(shift, index=oracle.adata.obs_names).groupby(
        oracle.adata.obs['cell_type'].values, observed=True).mean().sort_values(ascending=False)

    print(f'Mean shift magnitude by cell type after {perturbed_tf} KO:')
    print(mean_shift)

    if expected_affected_cluster:
        if expected_affected_cluster in mean_shift.index[:3]:
            print(f'OK: {expected_affected_cluster} among top affected clusters')
        else:
            print(f'WARNING: {expected_affected_cluster} not among top affected')

    return mean_shift
```

## Complete Pipeline Script

```python
import scanpy as sc
import pandas as pd
from arboreto.algo import grnboost2
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
from pyscenic.utils import modules_from_adjacencies
from ctxcore.rnkdb import FeatherRankingDatabase

def run_scenic_pipeline(adata_path, tf_list_path, db_paths, motif_annotations_path, output_prefix):
    '''Run complete pySCENIC pipeline.'''
    adata = sc.read_h5ad(adata_path)

    expr_matrix = pd.DataFrame(
        adata.raw.X.toarray() if hasattr(adata.raw.X, 'toarray') else adata.raw.X,
        index=adata.obs_names, columns=adata.raw.var_names
    )

    tf_names = pd.read_csv(tf_list_path, header=None)[0].tolist()
    tf_names = [tf for tf in tf_names if tf in expr_matrix.columns]

    print(f'Step 1: GRN inference with {len(tf_names)} TFs')
    adjacencies = grnboost2(expr_matrix, tf_names=tf_names, seed=42, verbose=True)

    print('Step 2: Regulon pruning')
    dbs = [FeatherRankingDatabase(db, name=os.path.splitext(os.path.basename(db))[0]) for db in db_paths]
    modules = list(modules_from_adjacencies(adjacencies, expr_matrix))
    df_motifs = prune2df(dbs, modules, motif_annotations_path, rank_threshold=5000, num_workers=8)
    regulons = df2regulons(df_motifs)
    print(f'Discovered {len(regulons)} regulons')

    print('Step 3: AUCell scoring')
    auc_matrix = aucell(expr_matrix, regulons, num_workers=8)
    adata.obsm['X_aucell'] = auc_matrix.loc[adata.obs_names].values
    adata.uns['regulon_names'] = [r.name for r in regulons]

    adata.write(f'{output_prefix}_scenic.h5ad')
    auc_matrix.to_csv(f'{output_prefix}_aucell.csv')

    print(f'Pipeline complete: {len(regulons)} regulons, AUCell matrix saved')
    return adata, regulons, auc_matrix
```

## Parameter Recommendations

| Step | Parameter | Recommendation |
|------|-----------|----------------|
| GRNBoost2 | min_targets | 10 (minimum targets per TF module) |
| RcisTarget | NES threshold | 3.0 (standard), 2.5 (permissive) |
| RcisTarget | databases | Use both 500bp and 10kbp upstream databases |
| AUCell | auc_threshold | 0.05 (fraction of ranked genes) |
| cisTopic | n_topics | Test 2x expected cell types |
| CellOracle | n_propagation | 3 (default signal propagation steps) |
| CellOracle | k (imputation) | int(0.025 * n_cells) (CellOracle tutorial rule; ~1250 at 50k cells) |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Near-empty regulons | Species/namespace/DB-vintage mismatch across TF-list, ranking DB, motif2TF | Pin all three to the SAME species + assembly + collection vintage |
| "Hybrid-state regulator" artifact | Ran GRN on a doublet-contaminated or un-integrated object | Infer on cleaned, doublet-free, batch-controlled cells; run SCENIC once on the integrated object |
| Inflated adjacencies / everything correlates | Inferred on imputed/smoothed counts | Use RAW counts; imputation only inside CellOracle's simulation scope |
| Regulon "validated" by TF-expression correlation | AUCell activity <-> TF mRNA circularity | Validate with an orthogonal modality (perturbation/ChIP), not the TF's own mRNA |
| ctx step returns empty | Missing/mismatched motif2TF annotation (most common) | Confirm the `.tbl` matches the DB vintage + species |
| SCENIC+ peaks miss rare types | Called peaks before/without cell-type labels | Label cells first; pycisTopic calls per-celltype pseudobulk peaks |
| Perturbation magnitudes reported as quantitative | Over-read the direction-only local model | Report direction + a baseline; never a quantitative KO magnitude |
| Modules called "regulons" without directionality | Skipped the cisTarget ctx pruning step | Run ctx; co-expression modules become regulons only after motif pruning |
| < 50 regulons / > 500 regulons | Strict pruning-wrong TF list / permissive thresholds | Lower NES to 2.5 (verify species) / raise NES to 3.5 |
| GRNBoost2 hangs or memory error | arboreto dask backend / large dataset | Use bundled multiprocessing; subsample to ~50k cells for GRNBoost2 |

## References

- Van de Sande B, Flerin C, Davie K, et al (2020) A scalable SCENIC workflow for single-cell gene regulatory network analysis. *Nature Protocols* 15:2247-2276. DOI 10.1038/s41596-020-0336-2. (pySCENIC 3-step; multi-run stability.)
- Bravo González-Blas C, De Winter S, Hulselmans G, et al (2023) SCENIC+: single-cell multiomic inference of enhancers and gene regulatory networks. *Nature Methods* 20:1355-1367. DOI 10.1038/s41592-023-01938-4. (eRegulons; needs paired multiome + cell-type labels before peak calling.)
- Kamimoto K, Stringa B, Hoffmann CM, et al (2023) Dissecting cell identity via network inference and in silico gene perturbation. *Nature* 614:742-751. DOI 10.1038/s41586-022-05688-9. (CellOracle; direction-only in-silico perturbation.)

## Related Skills

- gene-regulatory-networks/scenic-regulons - pySCENIC implementation details
- gene-regulatory-networks/multiomics-grn - SCENIC+ enhancer-driven GRNs
- gene-regulatory-networks/perturbation-simulation - CellOracle details
- single-cell/clustering - Upstream cell type annotation
- single-cell/preprocessing - QC and normalization before GRN inference
- atac-seq/single-cell-atac - scATAC preprocessing for SCENIC+ Multiome input
- atac-seq/co-accessibility - Cicero / SCENIC+ cis-regulatory connections
- atac-seq/enhancer-gene-linking - ABC / ENCODE-rE2G enhancer-gene mapping
- atac-seq/motif-deviation - chromVAR for TF motif accessibility
- workflows/scrnaseq-pipeline - Upstream: provides the cleaned, annotated RNA object for Path A (pySCENIC)
- workflows/multiome-pipeline - Upstream: provides the paired RNA+ATAC object for Path B (SCENIC+)
