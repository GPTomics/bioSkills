# Reference: scanpy 1.10+ | Verify API if version differs
# Hashtag demultiplexing with hashsolo: Bayesian sample assignment that handles few hashes
# and many negatives. HTO counts must live as columns in adata.obs.

import pandas as pd
import scanpy as sc
import scanpy.external as sce

# adata: AnnData with the GEX matrix; hto_counts_df: per-cell raw HTO counts indexed by cell barcode
# (from the CITE-seq/HTO count matrix, e.g. CellRanger's Antibody Capture output written to CSV)
adata = sc.read_h5ad('gex.h5ad')
hto_counts_df = pd.read_csv('hto_counts.csv', index_col=0)
hto_cols = ['HTO_A', 'HTO_B', 'HTO_C', 'HTO_D']
adata.obs[hto_cols] = hto_counts_df.loc[adata.obs_names, hto_cols]

# priors are ordered [negative, singlet, doublet]; raise the doublet prior for higher 10x loading
sce.pp.hashsolo(adata, cell_hashing_columns=hto_cols, priors=(0.01, 0.8, 0.19))

classification = adata.obs['Classification'].value_counts()

# Cross-sample doublet fraction calibrates the expected total doublet rate; sanity-check it
doublet_rate = (adata.obs['Classification'] == 'Doublet').mean()

singlets = adata[~adata.obs['Classification'].isin(['Negative', 'Doublet'])].copy()

classification
doublet_rate
singlets
