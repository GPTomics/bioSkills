'''Vertical (same-pixel) spatial multi-omics integration demo: RNA + protein co-profiled
on the SAME pixels (spatial CITE-seq / DBiT-seq style), joined with MOFA.

Serial-section modalities are different cells -- this toy is deliberately SAME-PIXEL so a
joint same-cell method (MOFA) is valid. For adjacent sections, register instead (PASTE).
'''
# Reference: muon 0.1+, mofapy2 0.7+, scanpy 1.10+, anndata 0.10+ | Verify API if version differs
import numpy as np
import anndata as ad
import mudata as md
import scanpy as sc
import muon as mu

rng = np.random.default_rng(0)
n_pixels, n_genes, n_proteins, n_groups = 400, 200, 30, 3

# shared latent groups drive BOTH modalities -- this is the cross-modal structure MOFA should recover
group = rng.integers(0, n_groups, size=n_pixels)
gene_means = rng.uniform(0.5, 4.0, size=(n_groups, n_genes))
prot_means = rng.uniform(0.5, 4.0, size=(n_groups, n_proteins))

rna_counts = rng.poisson(gene_means[group]).astype('float32')
prot_counts = rng.poisson(prot_means[group]).astype('float32')
coords = rng.uniform(0, 1000, size=(n_pixels, 2))                  # one shared pixel coordinate frame

pixel_ids = np.array([f'px{i}' for i in range(n_pixels)])
rna = ad.AnnData(rna_counts, obs={'pixel': pixel_ids})
prot = ad.AnnData(prot_counts, obs={'pixel': pixel_ids})
rna.obs_names = pixel_ids
prot.obs_names = pixel_ids
rna.var_names = [f'gene{j}' for j in range(n_genes)]
prot.var_names = [f'ADT{j}' for j in range(n_proteins)]
rna.obsm['spatial'] = coords

mdata = md.MuData({'rna': rna, 'prot': prot})
mu.pp.intersect_obs(mdata)                                         # keep only pixels present in BOTH modalities

sc.pp.normalize_total(mdata['rna'])
sc.pp.log1p(mdata['rna'])
sc.pp.highly_variable_genes(mdata['rna'])
mu.prot.pp.clr(mdata['prot'])                                      # protein is a panel: CLR, not log1p-of-counts

mu.tl.mofa(mdata, n_factors=5, use_var='highly_variable', outfile=None, verbose=False)

sc.pp.neighbors(mdata, use_rep='X_mofa')
sc.tl.leiden(mdata, flavor='igraph', n_iterations=2, directed=False)

mdata.obsm['X_mofa'].shape
n_clusters = mdata.obs['leiden'].nunique()
print(f'joint embedding {mdata.obsm["X_mofa"].shape}, leiden clusters: {n_clusters}, true groups: {n_groups}')
