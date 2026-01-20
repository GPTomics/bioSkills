#!/usr/bin/env python3
'''Batch integration with Harmony in Scanpy'''

import scanpy as sc
import harmonypy as hm
import matplotlib.pyplot as plt
import sys

def integrate_with_harmony(h5ad_files, output_prefix='integrated'):
    '''Integrate multiple scRNA-seq datasets with Harmony'''
    print('Loading datasets...')
    adatas = []
    for f in h5ad_files:
        adata = sc.read_h5ad(f)
        adata.obs['batch'] = f.replace('.h5ad', '').split('/')[-1]
        adatas.append(adata)

    print('Concatenating...')
    adata = sc.concat(adatas, label='batch', keys=[a.obs['batch'].iloc[0] for a in adatas])

    print('Preprocessing...')
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='batch')
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50)

    print('Running Harmony...')
    ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'batch')
    adata.obsm['X_pca_harmony'] = ho.Z_corr.T

    print('Post-integration...')
    sc.pp.neighbors(adata, use_rep='X_pca_harmony')
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5)

    print('Plotting...')
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    sc.pl.umap(adata, color='batch', ax=axes[0], show=False, title='By Batch')
    sc.pl.umap(adata, color='leiden', ax=axes[1], show=False, title='By Cluster')
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_umap.png', dpi=150)

    adata.write(f'{output_prefix}.h5ad')
    print(f'Saved: {output_prefix}.h5ad')

    return adata

if __name__ == '__main__':
    if len(sys.argv) > 1:
        integrate_with_harmony(sys.argv[1:])
    else:
        print('Usage: python harmony_integration.py sample1.h5ad sample2.h5ad sample3.h5ad')
