'''Fork-aware QC and normalization for spatial transcriptomics.

Demonstrates the platform-class fork on Squidpy built-in datasets: a Visium spot
section (library size tracks cellularity) and a MERFISH imaging section (low-count
targeted panel). Run: python preprocess_spatial.py
'''
# Reference: squidpy 1.5+, scanpy 1.10+, anndata 0.10+ | Verify API if version differs

import numpy as np
import scanpy as sc
import squidpy as sq


def negative_control_fdr(adata):
    '''FDR = mean counts per control feature / mean counts per real gene (imaging specificity proxy).'''
    ctrl = adata.var['control'].values
    if not ctrl.any():
        return float('nan')
    mean_ctrl = np.asarray(adata[:, ctrl].X.sum(axis=0)).ravel().mean()
    mean_gene = np.asarray(adata[:, ~ctrl].X.sum(axis=0)).ravel().mean()
    return mean_ctrl / mean_gene if mean_gene else float('nan')


def preprocess_spot(adata):
    '''Visium branch: UMI/gene/mito floors are tissue-dependent (OSTA DLPFC values, not law).'''
    adata.var['mt'] = adata.var_names.str.startswith(('MT-', 'mt-'))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, inplace=True)
    print(f'  spot UMI median {adata.obs.total_counts.median():.0f}, genes median {adata.obs.n_genes_by_counts.median():.0f}')
    sc.pp.filter_cells(adata, min_counts=600)
    sc.pp.filter_genes(adata, min_cells=10)
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_total(adata)            # library size carries biology -- this can blur spatial domains
    sc.pp.log1p(adata)
    return adata


def preprocess_imaging(adata):
    '''Imaging branch: low transcript floor, control FDR, area-style normalization (not gene-count-based).'''
    ctrl_prefixes = ('NegControlProbe', 'NegControlCodeword', 'BLANK', 'Blank', 'NegPrb')
    adata.var['control'] = adata.var_names.str.startswith(ctrl_prefixes)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['control'], percent_top=None, inplace=True)
    print(f'  imaging transcripts/cell median {adata.obs.total_counts.median():.0f} (tens-low hundreds, not thousands)')
    print(f'  negative-control FDR {negative_control_fdr(adata):.4f} (band ~<=0.01-0.05)')
    adata = adata[:, ~adata.var['control'].values].copy()
    sc.pp.filter_cells(adata, min_counts=10)   # Squidpy/community convention, NOT a vendor spec
    sc.pp.filter_genes(adata, min_cells=5)
    adata.layers['counts'] = adata.X.copy()
    if 'volume' in adata.obs:                  # area/volume denominator is panel-composition-independent (Moffitt-style)
        sf = adata.obs['volume'].values / adata.obs['volume'].median()
        adata.X = adata.X / sf[:, None]
    sc.pp.log1p(adata)
    return adata


visium = sq.datasets.visium_hne_adata()
print(f'Visium spot section: {visium.n_obs} spots, {visium.n_vars} genes')
visium = preprocess_spot(visium)
print(f'  after: {visium.n_obs} spots, {visium.n_vars} genes\n')

merfish = sq.datasets.merfish()
print(f'MERFISH imaging section: {merfish.n_obs} cells, {merfish.n_vars} genes')
merfish = preprocess_imaging(merfish)
print(f'  after: {merfish.n_obs} cells, {merfish.n_vars} genes')
