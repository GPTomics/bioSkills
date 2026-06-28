'''Spatially variable genes are mostly cell-type markers: an SVG-vs-HVG demonstration.'''
# Reference: squidpy 1.4+, scanpy 1.10+, anndata 0.10+ | Verify API if version differs
# Runs against a built-in Squidpy dataset; writes no files.

import squidpy as sq
import scanpy as sc

adata = sq.datasets.visium_hne_adata()       # Visium H&E mouse brain with 'cluster' labels
print(f'Loaded {adata.n_obs} spots, {adata.n_vars} genes')

sc.pp.highly_variable_genes(adata, n_top_genes=2000)
hvg = set(adata.var_names[adata.var['highly_variable']])

# Test SVGs over a broader universe than HVG so the SVG-not-HVG subset can exist at all
sc.pp.filter_genes(adata, min_cells=50)
universe = adata.var_names.tolist()

# The graph IS the model: n_neighs=6 mimics the Visium hex lattice (see spatial-neighbors)
sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=6)

# Analytic normal p-value + FDR over the full universe (not just HVGs) so overlap is meaningful
sq.gr.spatial_autocorr(adata, mode='moran', genes=universe, corr_method='fdr_bh', n_jobs=1)
moran = adata.uns['moranI']

# Threshold on effect size AND FDR, not p alone -- large n makes trivial autocorrelation 'significant'
svg = moran[(moran['I'] > 0.1) & (moran['pval_norm_fdr_bh'] < 0.05)]
svg_genes = set(svg.index)
print(f'\n{len(svg_genes)} spatially variable genes (I > 0.1, FDR < 0.05)')

# The central trap: SVGs overlap heavily with HVGs because both track spatially clustered cell types
overlap = len(svg_genes & hvg)
print(f'{overlap} of {len(svg_genes)} SVGs ({100 * overlap / max(len(svg_genes), 1):.0f}%) are also HVGs')
print('The SVG-not-HVG subset is where spatial information actually lives:')
print(sorted(svg_genes - hvg)[:15])

print('\nTop SVGs (inspect: are these cell-type markers rather than spatially regulated genes?):')
print(svg.sort_values('I', ascending=False).head(10)[['I', 'pval_norm_fdr_bh']])
