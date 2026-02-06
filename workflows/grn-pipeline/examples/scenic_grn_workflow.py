'''pySCENIC RNA-only GRN inference pipeline.'''

import scanpy as sc
import pandas as pd
import numpy as np
from arboreto.algo import grnboost2
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
from ctxcore.rnkdb import FeatherRankingDatabase
from ctxcore.genesig import GeneSignature
import matplotlib.pyplot as plt

# Configuration
ADATA_PATH = 'processed.h5ad'
TF_LIST = 'allTFs_hg38.txt'
DB_PATHS = [
    'hg38_500bp_up_100bp_down.genes_vs_motifs.rankings.feather',
    'hg38_10kbp_up_10kbp_down.genes_vs_motifs.rankings.feather'
]
MOTIF_ANNOTATIONS = 'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl'
OUTPUT_PREFIX = 'scenic_results'
NUM_WORKERS = 8

# Load data
# Example: Use any processed scRNA-seq AnnData (e.g., from scanpy PBMC tutorial)
adata = sc.read_h5ad(ADATA_PATH)
print(f'Loaded: {adata.n_obs} cells, {adata.n_vars} genes')

# Extract raw count expression matrix (GRNBoost2 works on counts)
expr_matrix = pd.DataFrame(
    adata.raw.X.toarray() if hasattr(adata.raw.X, 'toarray') else adata.raw.X,
    index=adata.obs_names, columns=adata.raw.var_names
)

# Load TF list and filter to genes in dataset
# Human TFs: https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt
# Mouse TFs: https://resources.aertslab.org/cistarget/tf_lists/allTFs_mm.txt
tf_names = pd.read_csv(TF_LIST, header=None)[0].tolist()
tf_names = [tf for tf in tf_names if tf in expr_matrix.columns]
print(f'TFs in dataset: {len(tf_names)}')

# Step 1: GRN inference with GRNBoost2
# Subsample if > 50k cells to manage memory
if expr_matrix.shape[0] > 50000:
    sample_idx = np.random.choice(expr_matrix.shape[0], 50000, replace=False)
    expr_sub = expr_matrix.iloc[sample_idx]
    print(f'Subsampled to {len(sample_idx)} cells for GRN inference')
else:
    expr_sub = expr_matrix

adjacencies = grnboost2(expr_sub, tf_names=tf_names, seed=42, verbose=True)
adjacencies.to_csv(f'{OUTPUT_PREFIX}_adjacencies.tsv', sep='\t', index=False)
print(f'TF-target adjacencies: {len(adjacencies)}')

# Step 2: Regulon pruning with RcisTarget
# cisTarget databases: ~10 GB each, download once from:
# https://resources.aertslab.org/cistarget/databases/
dbs = [FeatherRankingDatabase(db) for db in DB_PATHS]

# Build modules from adjacencies (minimum 10 targets per TF)
modules = [GeneSignature(name=tf, gene2weight=dict(zip(grp['target'], grp['importance'])))
           for tf, grp in adjacencies.groupby('TF') if len(grp) >= 10]
print(f'TF modules before pruning: {len(modules)}')

# NES threshold 3.0: Standard motif enrichment cutoff
# Lower to 2.5 for more permissive; raise to 3.5 for stringent
df_motifs = prune2df(dbs, modules, MOTIF_ANNOTATIONS, num_workers=NUM_WORKERS)
regulons = df2regulons(df_motifs)
print(f'Regulons after pruning: {len(regulons)}')

# QC: Check regulon count in expected range (50-500)
if len(regulons) < 50:
    print('WARNING: Few regulons. Lower NES threshold or check TF list.')
elif len(regulons) > 500:
    print('WARNING: Many regulons. Consider stricter NES threshold.')

# Step 3: AUCell regulon activity scoring
# Uses full expression matrix (all cells)
auc_matrix = aucell(expr_matrix, regulons, num_workers=NUM_WORKERS)
adata.obsm['X_aucell'] = auc_matrix.loc[adata.obs_names].values
adata.uns['regulon_names'] = [r.name for r in regulons]

# Save results
auc_matrix.to_csv(f'{OUTPUT_PREFIX}_aucell.csv')
adata.write(f'{OUTPUT_PREFIX}_scenic.h5ad')

# QC: Check known lineage TFs
known_tfs = ['PAX6', 'SOX2', 'GATA1', 'SPI1', 'FOXP3', 'TBX21', 'EBF1', 'IRF4']
regulon_names = [r.name for r in regulons]
found_known = [tf for tf in known_tfs if tf in regulon_names]
print(f'Known lineage TFs found: {found_known}')

# Visualization: top regulons by specificity
if 'cell_type' in adata.obs.columns:
    from scipy.stats import mannwhitneyu
    cell_types = adata.obs['cell_type'].unique()

    regulon_specificity = {}
    for i, rname in enumerate(regulon_names):
        scores = auc_matrix.iloc[:, i]
        max_diff = 0
        for ct in cell_types:
            ct_mask = adata.obs['cell_type'] == ct
            other_mask = ~ct_mask
            stat, pval = mannwhitneyu(scores[ct_mask], scores[other_mask], alternative='greater')
            if pval < 0.01:
                diff = scores[ct_mask].mean() - scores[other_mask].mean()
                max_diff = max(max_diff, diff)
        regulon_specificity[rname] = max_diff

    top_regulons = sorted(regulon_specificity, key=regulon_specificity.get, reverse=True)[:20]
    print(f'Top 20 cell-type-specific regulons: {top_regulons}')

print(f'\nPipeline complete: {len(regulons)} regulons discovered')
print(f'AUCell matrix: {OUTPUT_PREFIX}_aucell.csv')
print(f'Annotated AnnData: {OUTPUT_PREFIX}_scenic.h5ad')
