'''SCENIC+ multiome GRN inference pipeline.'''
# Reference: anndata 0.10+, pandas 2.2+, scanpy 1.10+, scipy 1.12+ | Verify API if version differs

import scanpy as sc
import numpy as np
from pycisTopic.cistopic_class import create_cistopic_object
from pycisTopic.lda_models import run_cgs_models, evaluate_models
from pycisTopic.topic_binarization import binarize_topics
from pycistarget.motif_enrichment_cistarget import run_cistarget
from scenicplus.scenicplus_class import create_SCENICPLUS_object
from scenicplus.grn_builder.gsea_approach import build_grn
import celloracle as co

# Configuration
RNA_PATH = 'rna_processed.h5ad'
ATAC_PATH = 'atac_processed.h5ad'
OUTPUT_PREFIX = 'scenicplus_results'
NUM_WORKERS = 8

# Load multiome data (paired RNA + ATAC)
# Example: 10x Multiome dataset
adata_rna = sc.read_h5ad(RNA_PATH)
adata_atac = sc.read_h5ad(ATAC_PATH)
print(f'RNA: {adata_rna.n_obs} cells, {adata_rna.n_vars} genes')
print(f'ATAC: {adata_atac.n_obs} cells, {adata_atac.n_vars} regions')

# Step 1: ATAC topic modeling with cisTopic
cistopic_obj = create_cistopic_object(
    fragment_matrix=adata_atac.X,
    cell_names=adata_atac.obs_names.tolist(),
    region_names=adata_atac.var_names.tolist()
)

# Test multiple topic numbers around 2x expected cell types
# n_iter=300: Standard for convergence; increase to 500 for complex datasets
n_clusters = len(adata_rna.obs['cell_type'].unique()) if 'cell_type' in adata_rna.obs else 10
topic_range = [n_clusters, n_clusters * 2, n_clusters * 3, n_clusters * 4]
print(f'Testing topic counts: {topic_range}')

models = run_cgs_models(
    cistopic_obj,
    n_topics=topic_range,
    n_cpu=NUM_WORKERS,
    n_iter=300,
    random_state=42
)

# Select best model by log-likelihood
model = evaluate_models(models, select_model=True)
cistopic_obj.add_LDA_model(model)
print(f'Selected model with {model.n_topic} topics')

# Step 2: Binarize topics to get region sets
region_bin = binarize_topics(cistopic_obj, method='otsu')
region_sets = {topic: regions for topic, regions in region_bin.items()}
print(f'Binarized {len(region_sets)} topics into region sets')

# Step 3: Motif enrichment on accessible regions
# NES threshold 3.0: Standard; auc_threshold 0.005: Top 0.5% of regions
cistarget_results = run_cistarget(
    region_sets=region_sets,
    species='homo_sapiens',
    auc_threshold=0.005,
    nes_threshold=3.0,
    rank_threshold=0.05,
    n_cpu=NUM_WORKERS
)

# Step 4: Build enhancer-driven GRN with SCENIC+
scplus_obj = create_SCENICPLUS_object(
    adata_rna=adata_rna,
    cistopic_obj=cistopic_obj,
    menr=cistarget_results
)

# min_target_genes=10: Minimum targets per eRegulon
# adj_pval_thr=0.1: FDR threshold for region-gene links
build_grn(
    scplus_obj,
    min_target_genes=10,
    adj_pval_thr=0.1,
    min_regions_per_gene=0
)

n_regulons = len(scplus_obj.uns['eRegulon_AUC'].columns)
print(f'eGRNs discovered: {n_regulons} enhancer-driven regulons')

# QC: Check regulon count
if n_regulons < 50:
    print('WARNING: Few eRegulons. Lower NES threshold or adj_pval_thr.')
elif n_regulons > 500:
    print('WARNING: Many eRegulons. Consider stricter thresholds.')

# Step 5 (optional): CellOracle perturbation simulation
oracle = co.Oracle()
oracle.import_anndata_as_raw_count(
    adata=adata_rna,
    cluster_column_name='cell_type',
    embedding_name='X_umap'
)

# Import GRN links from SCENIC+ adjacencies
links = co.utility.load_links(scplus_obj.uns.get('eRegulon_metadata', {}))
oracle.addTFinfo_dictionary(links)
oracle.perform_PCA()
# k=30: Neighborhood size for imputation; increase for smoother results
oracle.knn_imputation(k=30)

# Simulate TF knockout example
# Perturb MYC to 0 (complete knockout)
# n_propagation=3: Signal propagation steps through GRN
oracle.simulate_shift(perturb_condition={'MYC': 0.0}, n_propagation=3)
oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
oracle.calculate_embedding_shift(sigma_corr=0.05)

# Report perturbation scores by cell type
ps = oracle.perturbation_scores
if 'cell_type' in oracle.adata.obs:
    mean_shift = ps.groupby(oracle.adata.obs['cell_type']).mean()
    print('Perturbation scores by cell type (MYC KO):')
    print(mean_shift.sort_values(ascending=False))

print(f'\nPipeline complete: {n_regulons} enhancer-driven regulons')
print(f'Results stored in SCENIC+ object')
