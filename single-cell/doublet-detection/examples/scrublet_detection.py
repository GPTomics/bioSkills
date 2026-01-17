import scrublet as scr
import scanpy as sc
import matplotlib.pyplot as plt
import sys

input_path = sys.argv[1] if len(sys.argv) > 1 else 'filtered_feature_bc_matrix/'

adata = sc.read_10x_mtx(input_path)
print(f'Loaded {adata.n_obs} cells')

scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.06)
doublet_scores, predicted_doublets = scrub.scrub_doublets()

adata.obs['doublet_score'] = doublet_scores
adata.obs['is_doublet'] = predicted_doublets

n_doublets = predicted_doublets.sum()
pct_doublets = 100 * predicted_doublets.mean()
print(f'Detected {n_doublets} doublets ({pct_doublets:.1f}%)')

scrub.plot_histogram()
plt.savefig('doublet_histogram.pdf')
plt.close()

adata_clean = adata[~adata.obs['is_doublet']].copy()
print(f'Kept {adata_clean.n_obs} singlets')

adata_clean.write('singlets.h5ad')
