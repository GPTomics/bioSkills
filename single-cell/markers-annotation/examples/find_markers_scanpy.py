'''Find marker genes with Scanpy'''

import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('clustered.h5ad')

sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')

sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False, save='_markers.png')

markers = sc.get.rank_genes_groups_df(adata, group=None)
print('Top markers per cluster:')
print(markers.groupby('group').head(5)[['group', 'names', 'logfoldchanges', 'pvals_adj']])

markers.to_csv('all_markers.csv', index=False)

pbmc_markers = ['CD3D', 'CD8A', 'MS4A1', 'CD14', 'FCGR3A', 'NKG7']
sc.pl.dotplot(adata, var_names=pbmc_markers, groupby='leiden', save='_dotplot.png')

cluster_annotations = {'0': 'T cells', '1': 'Monocytes', '2': 'B cells'}
adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_annotations).fillna('Unknown')
sc.pl.umap(adata, color='cell_type', save='_celltypes.png')

adata.write_h5ad('annotated.h5ad')
print('Saved annotated data')
