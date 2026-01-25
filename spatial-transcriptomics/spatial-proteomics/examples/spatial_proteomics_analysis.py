#!/usr/bin/env python3
"""
Spatial Proteomics Analysis with scimap
Demonstrates analysis of CODEX/IMC/MIBI data
"""

import scimap as sm
import anndata as ad
import pandas as pd
import numpy as np

# =============================================================================
# Load and Preprocess Data
# =============================================================================

def load_spatial_proteomics(filepath):
    """Load spatial proteomics data.

    Expected format: h5ad with
    - adata.X: cell x marker intensity matrix
    - adata.obsm['spatial']: x,y coordinates
    - adata.obs: cell metadata (optional: fov, region, etc.)
    """
    adata = ad.read_h5ad(filepath)

    # Ensure spatial coordinates exist
    if 'spatial' not in adata.obsm:
        if 'X_centroid' in adata.obs and 'Y_centroid' in adata.obs:
            adata.obsm['spatial'] = adata.obs[['X_centroid', 'Y_centroid']].values
        else:
            raise ValueError('Spatial coordinates not found')

    return adata


def preprocess(adata):
    """Standard preprocessing for spatial proteomics."""
    # Log transform
    # Stabilizes variance, reduces impact of outliers
    sm.pp.log1p(adata)

    # Rescale each marker to 0-1
    # Makes markers comparable for gating
    sm.pp.rescale(adata, gate=0.5)

    return adata


# =============================================================================
# Cell Phenotyping
# =============================================================================

def phenotype_cells_gating(adata):
    """Phenotype cells using marker gating.

    Gate values after rescale:
    - 0.5 is middle of distribution
    - Higher = more stringent
    """
    # Define phenotypes by marker combinations
    phenotype_markers = {
        'T_cell': ['CD3+', 'CD45+'],
        'T_helper': ['CD3+', 'CD4+', 'CD45+'],
        'T_cytotoxic': ['CD3+', 'CD8+', 'CD45+'],
        'B_cell': ['CD20+', 'CD45+'],
        'Macrophage': ['CD68+', 'CD163+'],
        'Dendritic': ['CD11c+', 'HLA-DR+'],
        'Tumor': ['panCK+'],
        'Proliferating': ['Ki67+']
    }

    # Apply phenotyping
    # gate=0.5 means marker must be above 0.5 (after rescale)
    sm.tl.phenotype_cells(
        adata,
        phenotype=phenotype_markers,
        gate=0.5,
        label='phenotype'
    )

    print('Phenotype distribution:')
    print(adata.obs['phenotype'].value_counts())

    return adata


def phenotype_cells_clustering(adata, resolution=1.0):
    """Phenotype cells using unsupervised clustering."""
    # Leiden clustering
    sm.tl.cluster(adata, method='leiden', resolution=resolution)

    print('Cluster distribution:')
    print(adata.obs['leiden'].value_counts())

    # Annotate clusters based on marker expression
    # This requires manual inspection of cluster markers
    return adata


# =============================================================================
# Spatial Analysis
# =============================================================================

def spatial_neighbors(adata, method='knn', k=10):
    """Build spatial neighbor graph.

    Args:
        method: 'knn' or 'radius'
        k: number of neighbors (for knn)
    """
    sm.tl.spatial_distance(
        adata,
        x_coordinate='X_centroid',
        y_coordinate='Y_centroid'
    )
    return adata


def spatial_interaction_analysis(adata, phenotype_col='phenotype'):
    """Analyze cell-cell spatial interactions.

    Returns enrichment/depletion of cell type pairs
    compared to random distribution.
    """
    sm.tl.spatial_interaction(
        adata,
        phenotype=phenotype_col,
        method='knn',
        knn=10,
        permutation=1000  # For p-value calculation
    )

    # Results stored in adata.uns['spatial_interaction']
    return adata


def find_cellular_neighborhoods(adata, phenotype_col='phenotype'):
    """Identify recurring cellular neighborhoods."""
    # Spatial clustering based on cell type composition
    sm.tl.spatial_cluster(
        adata,
        phenotype=phenotype_col,
        method='kmeans',
        k=5  # Number of neighborhood types
    )

    print('Neighborhood distribution:')
    print(adata.obs['spatial_cluster'].value_counts())

    return adata


# =============================================================================
# Visualization
# =============================================================================

def visualize_results(adata):
    """Generate standard visualizations."""
    # Spatial scatter plot by phenotype
    sm.pl.spatial_scatterPlot(
        adata,
        colorBy='phenotype',
        x='X_centroid',
        y='Y_centroid',
        s=3,
        figsize=(10, 10)
    )

    # Spatial interaction heatmap
    if 'spatial_interaction' in adata.uns:
        sm.pl.spatial_interaction(adata, summarize_plot=True)

    # Neighborhood composition
    if 'spatial_cluster' in adata.obs:
        sm.pl.spatial_scatterPlot(
            adata,
            colorBy='spatial_cluster',
            x='X_centroid',
            y='Y_centroid'
        )


# =============================================================================
# Main Workflow
# =============================================================================

if __name__ == '__main__':
    # Example workflow
    print('=== Spatial Proteomics Analysis ===\n')

    # For demonstration, create synthetic data
    # In practice, load real CODEX/IMC data
    n_cells = 5000
    n_markers = 30

    # Create synthetic AnnData
    np.random.seed(42)
    adata = ad.AnnData(
        X=np.random.lognormal(0, 1, (n_cells, n_markers)),
        obs=pd.DataFrame({
            'X_centroid': np.random.uniform(0, 1000, n_cells),
            'Y_centroid': np.random.uniform(0, 1000, n_cells),
            'fov': np.random.choice(['FOV1', 'FOV2', 'FOV3'], n_cells)
        }),
        var=pd.DataFrame(index=[f'Marker_{i}' for i in range(n_markers)])
    )
    adata.obsm['spatial'] = adata.obs[['X_centroid', 'Y_centroid']].values

    print('1. Preprocessing...')
    adata = preprocess(adata)

    print('\n2. Clustering cells...')
    adata = phenotype_cells_clustering(adata, resolution=0.5)

    print('\n3. Building spatial graph...')
    adata = spatial_neighbors(adata)

    print('\n4. Analyzing spatial interactions...')
    adata = spatial_interaction_analysis(adata, phenotype_col='leiden')

    print('\n5. Finding cellular neighborhoods...')
    adata = find_cellular_neighborhoods(adata, phenotype_col='leiden')

    print('\n=== Analysis Complete ===')
    print(f'Total cells: {adata.n_obs}')
    print(f'Markers: {adata.n_vars}')
