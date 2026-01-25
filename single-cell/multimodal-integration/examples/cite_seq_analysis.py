import muon as mu
import scanpy as sc
import anndata as ad
import numpy as np

# Read 10X multiome data (RNA + ADT)
# muon reads 10X filtered_feature_bc_matrix folder with both modalities
mdata = mu.read_10x_h5('filtered_feature_bc_matrix.h5')

# Access individual modalities
rna = mdata.mod['rna']
adt = mdata.mod['prot']  # ADT/protein

# RNA QC
rna.var['mt'] = rna.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(rna, qc_vars=['mt'], percent_top=None, inplace=True)

# Filter cells - keep cells present in both modalities
mu.pp.intersect_obs(mdata)

# Filter by QC metrics
# nFeature_RNA > 200: minimum gene complexity
# nFeature_RNA < 5000: remove potential doublets
# percent.mt < 20: remove dying cells
rna_mask = (rna.obs['n_genes_by_counts'] > 200) & (rna.obs['n_genes_by_counts'] < 5000) & (rna.obs['pct_counts_mt'] < 20)
mdata = mdata[rna_mask].copy()

# RNA preprocessing
sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)
sc.pp.highly_variable_genes(rna, n_top_genes=2000)
sc.pp.scale(rna, max_value=10)
sc.tl.pca(rna, n_comps=30)

# ADT preprocessing (CLR normalization)
# margin=0: normalize across cells (standard for ADT)
mu.prot.pp.clr(adt, axis=0)
sc.pp.scale(adt, max_value=10)
# Use all ADT features for PCA (typically 10-200 markers)
sc.tl.pca(adt, n_comps=min(18, adt.n_vars - 1))

# Weighted nearest neighbors (WNN) for multimodal integration
mu.pp.neighbors(mdata, key_added='wnn')

# Clustering on WNN graph
sc.tl.leiden(mdata, resolution=0.5, key_added='wnn_clusters', neighbors_key='wnn')

# UMAP on WNN
sc.tl.umap(mdata, neighbors_key='wnn')

# Save UMAP plot
sc.settings.figdir = './'
sc.pl.umap(mdata, color='wnn_clusters', save='_cite_seq_wnn.pdf')

# Save processed data
mdata.write('cite_seq_analyzed.h5mu')
print('Analysis complete. Saved to cite_seq_analyzed.h5mu')
