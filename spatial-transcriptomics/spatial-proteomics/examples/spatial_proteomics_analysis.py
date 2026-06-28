#!/usr/bin/env python3
"""Spatial proteomics analysis: intensity transform + neighborhood enrichment.

Runs on squidpy's built-in IMC dataset (Jackson et al.), a real
multiplexed-imaging cell-by-marker matrix with continuous metal-channel
intensities, spatial coordinates, and a 'cell type' annotation. Demonstrates
the central reframe: intensity is continuous and confounded (arcsinh, NOT
log1p-of-counts), and spatial neighborhood results sit downstream of phenotyping.
"""
# Reference: squidpy 1.4+, scanpy 1.10+, anndata 0.10+ | Verify API if version differs

import numpy as np
import scanpy as sc
import squidpy as sq


def load_imc():
    """Load the built-in IMC dataset (~4668 cells x 34 metal channels)."""
    adata = sq.datasets.imc()
    # Intensities live in adata.X; 'cell type' holds the published phenotypes.
    return adata


def transform_intensities(adata, cofactor=5):
    """Variance-stabilize metal-channel intensities with arcsinh.

    cofactor=5 is the CyTOF convention, not auto-optimal for imaging -- too
    small over-expands near-zero noise into spurious positive populations, so
    it should be sanity-checked per dataset. This is NOT log1p-of-counts:
    intensity has no Poisson/NB count process.
    """
    adata.layers['intensity'] = adata.X.copy()
    adata.X = np.arcsinh(adata.X / cofactor)
    return adata


def neighborhood_enrichment(adata, cluster_key='cell type'):
    """Test which phenotypes are spatial neighbors more/less than chance.

    Imaging cells are a point cloud (not a grid) -> coord_type='generic'.
    Enrichment uses a label-permutation null; results inherit every
    phenotyping and segmentation error upstream.
    """
    sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=10)
    sq.gr.nhood_enrichment(adata, cluster_key=cluster_key)
    return adata


if __name__ == '__main__':
    adata = load_imc()
    print(f'Loaded IMC: {adata.n_obs} cells x {adata.n_vars} markers')
    print('Cell types:')
    print(adata.obs['cell type'].value_counts())

    adata = transform_intensities(adata, cofactor=5)
    print(f'\nIntensity range after arcsinh: {adata.X.min():.2f} to {adata.X.max():.2f}')

    adata = neighborhood_enrichment(adata)
    zscores = adata.uns['cell type_nhood_enrichment']['zscore']
    print(f'\nNeighborhood-enrichment z-score matrix: {zscores.shape}')
    print(f'Strongest co-localization z-score: {np.nanmax(zscores[~np.eye(len(zscores), dtype=bool)]):.1f}')
    print('\nAnalysis complete -- results are downstream of phenotyping and segmentation.')
