# Reference: scanpy 1.10+ | Verify API if version differs
import scanpy as sc
import matplotlib.pyplot as plt
import sys

input_path = sys.argv[1] if len(sys.argv) > 1 else 'filtered_feature_bc_matrix/'

adata = sc.read_10x_mtx(input_path)
print(f'Loaded {adata.n_obs} cells')

# Set the expected rate from recovered cells (~0.008 per 1000), not the flat 0.05 placeholder
expected_rate = 0.008 * adata.n_obs / 1000
# Maintained Scrublet path; for pooled data pass batch_key to run per sample, never on an integrated object
sc.pp.scrublet(adata, expected_doublet_rate=expected_rate)

n_doublets = int(adata.obs['predicted_doublet'].sum())
pct_doublets = 100 * adata.obs['predicted_doublet'].mean()
print(f'Detected {n_doublets} doublets ({pct_doublets:.1f}%)')

# Inspect the histogram: the auto-threshold needs a bimodal distribution; if unimodal, set a manual cutoff
sc.pl.scrublet_score_distribution(adata, show=False)
plt.savefig('doublet_histogram.pdf')
plt.close()

adata_clean = adata[~adata.obs['predicted_doublet']].copy()
print(f'Kept {adata_clean.n_obs} singlets')

adata_clean.write('singlets.h5ad')
