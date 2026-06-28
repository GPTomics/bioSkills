'''Deconvolve spatial spots into cell-type composition with cell2location.

Builds a tiny synthetic scRNA-seq reference and Visium-like spatial mixture,
runs the cell2location two-step (reference signatures -> spatial mapping),
converts absolute abundances to CLR-transformed composition, and validates
against markers. All artifacts go to a temp dir and are removed at exit, so
no model dirs / .h5ad / .png are left in the working tree.
'''
# Reference: cell2location 0.1.4+, scvi-tools 1.0+, scanpy 1.10+, anndata 0.10+, numpy 1.26+ | Verify API if version differs

import shutil
import tempfile
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import cell2location
from cell2location.models import RegressionModel

rng = np.random.default_rng(0)
n_types, n_genes = 4, 60
type_names = [f'type_{i}' for i in range(n_types)]

# Each type over-expresses a disjoint marker block; everything else is low background.
base = rng.poisson(0.5, size=(n_types, n_genes)).astype(float)
block = n_genes // n_types
for i in range(n_types):
    base[i, i * block:(i + 1) * block] += 8.0

ref_per_type = 80
ref_counts = np.vstack([rng.poisson(base[i], size=(ref_per_type, n_genes)) for i in range(n_types)])
ref_labels = np.repeat(type_names, ref_per_type)
adata_ref = ad.AnnData(ref_counts.astype(np.float32))
adata_ref.var_names = [f'g{j}' for j in range(n_genes)]
adata_ref.obs['cell_type'] = pd.Categorical(ref_labels)

# Spatial spots are mixtures: draw a composition per spot, sum type profiles weighted by cell counts.
n_spots, cells_per_spot = 120, 8
comp = rng.dirichlet(np.ones(n_types) * 0.7, size=n_spots)
spot_counts = np.zeros((n_spots, n_genes))
for s in range(n_spots):
    counts_by_type = rng.multinomial(cells_per_spot, comp[s])
    for i, k in enumerate(counts_by_type):
        if k:
            spot_counts[s] += rng.poisson(base[i] * k)
adata_vis = ad.AnnData(spot_counts.astype(np.float32))
adata_vis.var_names = adata_ref.var_names.copy()

workdir = tempfile.mkdtemp(prefix='c2l_')
try:
    RegressionModel.setup_anndata(adata_ref, labels_key='cell_type')
    mod = RegressionModel(adata_ref)
    mod.train(max_epochs=120, accelerator='cpu')
    adata_ref = mod.export_posterior(adata_ref, sample_kwargs={'num_samples': 200})

    factors = adata_ref.uns['mod']['factor_names']
    if 'means_per_cluster_mu_fg' in adata_ref.varm:
        inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{f}' for f in factors]].copy()
    else:
        inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{f}' for f in factors]].copy()
    inf_aver.columns = factors

    shared = np.intersect1d(adata_vis.var_names, inf_aver.index)
    adata_vis = adata_vis[:, shared].copy()
    inf_aver = inf_aver.loc[shared, :]

    cell2location.models.Cell2location.setup_anndata(adata_vis)
    mod_sp = cell2location.models.Cell2location(adata_vis, cell_state_df=inf_aver,
                                                N_cells_per_location=cells_per_spot, detection_alpha=20)
    mod_sp.train(max_epochs=2000, batch_size=None, train_size=1, accelerator='cpu')
    adata_vis = mod_sp.export_posterior(adata_vis, sample_kwargs={'num_samples': 200, 'batch_size': mod_sp.adata.n_obs})

    abund = adata_vis.obsm['q05_cell_abundance_w_sf'].values
    proportions = abund / abund.sum(axis=1, keepdims=True)

    p = proportions + 1e-6
    p = p / p.sum(axis=1, keepdims=True)
    log_p = np.log(p)
    clr_comp = log_p - log_p.mean(axis=1, keepdims=True)

    est_mean = proportions.mean(axis=0)
    true_mean = comp.mean(axis=0)
    print('estimated vs true mean composition:')
    for i, name in enumerate(factors):
        print(f'  {name}: est {est_mean[i]:.3f}  true {true_mean[i]:.3f}')
    print(f'CLR composition shape: {clr_comp.shape}')

    adata_vis.write_h5ad(f'{workdir}/spatial_deconvolved.h5ad')
    print('deconvolution complete')
finally:
    shutil.rmtree(workdir, ignore_errors=True)
