import scvelo as scv
import scanpy as sc
import numpy as np

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo')

adata = sc.read_h5ad('adata_clustered.h5ad')
ldata = scv.read('velocyto_output.loom')
adata = scv.utils.merge(adata, ldata)

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

scv.tl.recover_dynamics(adata, n_jobs=8)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)

scv.pl.velocity_embedding_stream(adata, basis='umap', color='clusters',
                                  save='velocity_stream.png', dpi=150)

scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', cmap='gnuplot', save='latent_time.png')

scv.tl.velocity_confidence(adata)
scv.pl.scatter(adata, color=['velocity_confidence', 'velocity_length'],
               save='velocity_confidence.png')

scv.tl.rank_velocity_genes(adata, groupby='clusters', min_corr=0.3)
velocity_genes = adata.uns['rank_velocity_genes']['names'][:10]
scv.pl.velocity(adata, var_names=list(velocity_genes.flatten()[:6]),
                basis='umap', save='top_velocity_genes.png')

adata.write('adata_with_velocity.h5ad')
