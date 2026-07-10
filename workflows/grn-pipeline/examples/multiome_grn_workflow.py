'''SCENIC+ multiome GRN: cisTopic + pycistarget prep, then the SCENIC+ Snakemake pipeline;
optional CellOracle perturbation off a motif-scan base GRN (NOT the SCENIC+ eRegulons).'''
# Reference: pycisTopic 2.0+, pycistarget 1.0+, SCENIC+ 1.0a1, CellOracle 0.18+, scanpy 1.10+ | Verify API if version differs

import glob
import subprocess
import scanpy as sc
import numpy as np
import pandas as pd
import pyranges as pr
from pycisTopic.cistopic_class import create_cistopic_object
from pycisTopic.lda_models import run_cgs_models, evaluate_models
from pycisTopic.topic_binarization import binarize_topics
from pycistarget.utils import region_names_to_coordinates
from pycistarget.motif_enrichment_cistarget import run_cistarget

# Configuration
RNA_PATH = 'rna_processed.h5ad'
ATAC_PATH = 'atac_processed.h5ad'
SCENICPLUS_DIR = 'scenicplus_run'
NUM_WORKERS = 8
# cisTarget ranking DB for the matching genome/region set; run_cistarget accepts the feather path
# directly and instantiates cisTargetDatabase itself. Prebuilt DBs: https://resources.aertslab.org/cistarget/
CTX_DB = '/path/to/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather'

# Load paired multiome data (RNA + ATAC in the SAME nuclei)
adata_rna = sc.read_h5ad(RNA_PATH)
adata_atac = sc.read_h5ad(ATAC_PATH)
print(f'RNA: {adata_rna.n_obs} cells, {adata_rna.n_vars} genes')
print(f'ATAC: {adata_atac.n_obs} cells, {adata_atac.n_vars} regions')

# Step 1: ATAC topic modeling with cisTopic (LDA)
cistopic_obj = create_cistopic_object(
    fragment_matrix=adata_atac.X.T,   # cisTopic wants regions x cells; AnnData .X is cells x regions
    cell_names=adata_atac.obs_names.tolist(),
    region_names=adata_atac.var_names.tolist()
)
# Test topic counts around 2x the number of cell types
n_clusters = len(adata_rna.obs['cell_type'].unique()) if 'cell_type' in adata_rna.obs else 10
topic_range = [n_clusters, n_clusters * 2, n_clusters * 3, n_clusters * 4]
models = run_cgs_models(cistopic_obj, n_topics=topic_range, n_cpu=NUM_WORKERS, n_iter=300, random_state=42)
model = evaluate_models(models, select_model=topic_range[2], return_model=True)   # int topic COUNT at the metric elbow, not True (True==1 -> non-existent 1-topic model)
cistopic_obj.add_LDA_model(model)

# Step 2: Binarize topics into region sets, then convert to PyRanges (what run_cistarget expects)
region_bin = binarize_topics(cistopic_obj, method='otsu')
region_sets = {topic: pr.PyRanges(region_names_to_coordinates(region_bin[topic].index.tolist()))
               for topic in region_bin}
print(f'Binarized {len(region_sets)} topics into region sets')

# Step 3: Motif enrichment on the accessible region sets
cistarget_results = run_cistarget(
    CTX_DB,                              # cisTarget ranking DB, required first positional arg
    region_sets=region_sets, specie='homo_sapiens',   # param is `specie`, not `species`
    auc_threshold=0.005, nes_threshold=3.0, rank_threshold=0.05, n_cpu=NUM_WORKERS
)

# Step 4: eGRN construction via the SCENIC+ SNAKEMAKE pipeline.
# The manual create_SCENICPLUS_object/build_grn object API is DEPRECATED in SCENIC+ 1.0 -- scaffold
# and run the Snakemake pipeline instead, pointing its config at the fragments, this RNA AnnData,
# the cell-type labels, and genome/motif-matched cisTarget DBs.
subprocess.run(['scenicplus', 'init_snakemake', '--out_dir', SCENICPLUS_DIR], check=True)
print(f'Edit {SCENICPLUS_DIR}/Snakemake/config/config.yaml, then run:')
print(f'  cd {SCENICPLUS_DIR}/Snakemake && snakemake --cores {NUM_WORKERS}')
# After the Snakemake run, read the direct (high-confidence) eRegulon table (filename is
# config-/version-dependent, so resolve it by glob):
matches = glob.glob(f'{SCENICPLUS_DIR}/**/eRegulon*direct*.tsv', recursive=True)
if matches:
    eregulons = pd.read_csv(matches[0], sep='\t')
    print(f'eRegulons: {eregulons["TF"].nunique()} enhancer-driven regulators')

# Step 5 (optional): CellOracle in-silico perturbation. The base GRN is a motif-scanned
# accessible-region scaffold (or a prebuilt CellOracle base GRN) -- NOT the SCENIC+ eRegulons.
# See gene-regulatory-networks/perturbation-simulation for base-GRN construction.
import celloracle as co
# A motif-scan scaffold over THESE ATAC peaks is preferable; the prebuilt promoter base GRN is the
# fallback when no matching scan exists. `version` must match the genome build of adata_rna.
base_grn = co.data.load_human_promoter_base_GRN()
oracle = co.Oracle()
oracle.import_anndata_as_raw_count(adata=adata_rna, cluster_column_name='cell_type', embedding_name='X_umap')
oracle.import_TF_data(TF_info_matrix=base_grn)
oracle.perform_PCA()
k = int(0.025 * oracle.adata.n_obs)
oracle.knn_imputation(n_pca_dims=50, k=k, balanced=True, b_sight=k * 8, b_maxl=k * 4)
links = oracle.get_links(cluster_name_for_GRN_unit='cell_type', alpha=10)
links.filter_links(p=0.001, weight='coef_abs', threshold_number=2000)
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)

oracle.simulate_shift(perturb_condition={'MYC': 0.0}, n_propagation=3)
oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
oracle.calculate_embedding_shift(sigma_corr=0.05)

# Shift magnitude per cell from the simulated embedding shift (direction, not calibrated magnitude)
shift = np.sqrt((oracle.delta_embedding ** 2).sum(axis=1))
mean_shift = pd.Series(shift, index=oracle.adata.obs_names).groupby(oracle.adata.obs['cell_type'].values, observed=True).mean().sort_values(ascending=False)   # observed=True: cell_type is categorical; the default keeps filtered-out categories as NaN rows
print('Mean shift magnitude by cell type after MYC KO:')
print(mean_shift)
