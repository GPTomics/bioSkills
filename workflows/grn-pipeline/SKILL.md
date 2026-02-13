---
name: bio-workflows-grn-pipeline
description: End-to-end gene regulatory network inference pipeline from processed single-cell data to regulon discovery and perturbation simulation. Supports RNA-only (pySCENIC) and multiome (SCENIC+) paths. Use when building gene regulatory networks from single-cell transcriptomic or multiome data.
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

Reference examples tested with: anndata 0.10+, pandas 2.2+, scanpy 1.10+, scipy 1.12+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Gene Regulatory Network Pipeline

**"Infer gene regulatory networks from my single-cell data"** → Orchestrate pySCENIC regulon inference (GRNBoost2, cisTarget, AUCell), CellOracle perturbation simulation, and regulon-based cell type characterization.

Complete workflow from processed single-cell data to regulon discovery and perturbation simulation.

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
from ctxcore.genesig import GeneSignature

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
from ctxcore.rnkdb import FeatherRankingDatabase

# cisTarget databases (~10 GB each, download once)
# Human: hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
# Mouse: mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
dbs = [FeatherRankingDatabase(db) for db in [
    'hg38_500bp_up_100bp_down.genes_vs_motifs.rankings.feather',
    'hg38_10kbp_up_10kbp_down.genes_vs_motifs.rankings.feather'
]]

motif_annotations = 'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl'

# Create modules from adjacencies (top 50 targets per TF)
modules = [GeneSignature(name=tf, gene2weight=dict(zip(grp['target'], grp['importance'])))
           for tf, grp in adjacencies.groupby('TF')
           if len(grp) >= 10]

# Prune modules using motif enrichment
# NES threshold 3.0: Standard; lower to 2.5 for more permissive
df_motifs = prune2df(dbs, modules, motif_annotations, num_workers=8)
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
    found = [tf for tf in known_tfs if tf in regulon_names]
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
    fragment_matrix=adata_atac.X,
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

# Select best model by log-likelihood
from pycisTopic.lda_models import evaluate_models
model = evaluate_models(models, select_model=True)
cistopic_obj.add_LDA_model(model)
```

### Step 2: Enhancer-TF Mapping

```python
from pycistarget.utils import region_names_to_coordinates
from pycistarget.motif_enrichment_cistarget import run_cistarget

region_sets = {}
from pycisTopic.topic_binarization import binarize_topics
region_bin = binarize_topics(cistopic_obj, method='otsu')

for topic in region_bin:
    region_sets[topic] = region_bin[topic]

# Run motif enrichment on accessible regions
cistarget_results = run_cistarget(
    region_sets=region_sets,
    species='homo_sapiens',
    auc_threshold=0.005,
    nes_threshold=3.0,
    rank_threshold=0.05,
    n_cpu=8
)
```

### Step 3: eGRN Construction

```python
from scenicplus.scenicplus_class import create_SCENICPLUS_object
from scenicplus.grn_builder.gsea_approach import build_grn

scplus_obj = create_SCENICPLUS_object(
    adata_rna=adata_rna,
    cistopic_obj=cistopic_obj,
    menr=cistarget_results
)

build_grn(
    scplus_obj,
    min_target_genes=10,
    adj_pval_thr=0.1,
    min_regions_per_gene=0
)

print(f'eGRNs: {len(scplus_obj.uns["eRegulon_AUC"].columns)} enhancer-driven regulons')
```

## CellOracle Perturbation Simulation

```python
import celloracle as co

oracle = co.Oracle()
oracle.import_anndata_as_raw_count(
    adata=adata,
    cluster_column_name='cell_type',
    embedding_name='X_umap'
)

# Import GRN (from pySCENIC adjacencies or SCENIC+ eGRNs)
links = co.utility.load_links(adjacencies)
oracle.addTFinfo_dictionary(links)

oracle.perform_PCA()
oracle.knn_imputation(k=30)

# Simulate TF knockout
oracle.simulate_shift(perturb_condition={'MYC': 0.0}, n_propagation=3)
oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)

# Perturbation score: magnitude of predicted shift
oracle.calculate_embedding_shift(sigma_corr=0.05)

# QC: shifts should match known biology (e.g., MYC KO reduces proliferation)
gradient = oracle.gradient_fitting(n_jobs=8)
```

### QC Checkpoint: Perturbation

```python
def validate_perturbation(oracle, perturbed_tf, expected_affected_cluster=None):
    '''
    QC gate: perturbation shifts should match known biology.
    - Transition probabilities should show directional shift
    - If expected_affected_cluster known, check it shows largest change
    '''
    ps = oracle.perturbation_scores
    mean_shift = ps.groupby(oracle.adata.obs['cell_type']).mean()

    print(f'Mean perturbation scores by cell type after {perturbed_tf} KO:')
    print(mean_shift.sort_values(ascending=False))

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
from ctxcore.rnkdb import FeatherRankingDatabase
from ctxcore.genesig import GeneSignature

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
    dbs = [FeatherRankingDatabase(db) for db in db_paths]
    modules = [GeneSignature(name=tf, gene2weight=dict(zip(grp['target'], grp['importance'])))
               for tf, grp in adjacencies.groupby('TF') if len(grp) >= 10]
    df_motifs = prune2df(dbs, modules, motif_annotations_path, num_workers=8)
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
| CellOracle | k (imputation) | 30 (neighborhood size for imputation) |

## Troubleshooting

| Issue | Likely Cause | Solution |
|-------|--------------|----------|
| < 50 regulons | Strict pruning or wrong TF list | Lower NES threshold to 2.5; verify TF list species |
| > 500 regulons | Permissive thresholds | Increase NES threshold to 3.5 |
| AUCell flat | Low regulon quality or normalization issue | Use raw counts; check regulon gene overlap |
| CellOracle no shift | Weak GRN or wrong TF | Verify TF expressed in target cells |
| Memory error | Large dataset | Subsample to 50k cells for GRNBoost2 |

## Related Skills

- gene-regulatory-networks/scenic-regulons - pySCENIC implementation details
- gene-regulatory-networks/multiomics-grn - SCENIC+ enhancer-driven GRNs
- gene-regulatory-networks/perturbation-simulation - CellOracle details
- single-cell/clustering - Upstream cell type annotation
- single-cell/preprocessing - QC and normalization before GRN inference
