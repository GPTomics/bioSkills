'''
CellOracle TF perturbation simulation workflow.

Requires:
- clustered.h5ad: preprocessed scRNA-seq with UMAP and cell type labels
- atac_peaks.bed: accessible chromatin regions (from scATAC or bulk ATAC)
'''

import scanpy as sc
import celloracle as co
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def build_base_grn(peaks_bed, ref_genome='hg38'):
    '''Scan accessible regions for TF motifs to build base GRN.'''
    peaks = pd.read_csv(peaks_bed, sep='\t', header=None, names=['chr', 'start', 'end'])

    tfi = co.motif_analysis.TFinfo(peak_data_frame=peaks, ref_genome=ref_genome)
    # fpr 0.02: false positive rate for motif scanning
    tfi.scan(fpr=0.02)
    tfi.filter_motifs_by_score(threshold=10)
    tfi.make_TFinfo_dataframe_and_target_gene_dataframe()

    base_grn = tfi.to_dataframe()
    base_grn.to_parquet('base_grn.parquet')
    print(f'Base GRN: {len(base_grn)} TF-target links, {base_grn["source"].nunique()} TFs')
    return base_grn


def construct_grn(adata, base_grn, cluster_column='cell_type'):
    '''Fit GRN from scRNA-seq expression and base GRN.'''
    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(
        adata=adata,
        cluster_column_name=cluster_column,
        embedding_name='X_umap'
    )
    oracle.import_TF_data(TF_info_matrix=base_grn)

    oracle.perform_PCA()
    oracle.knn_imputation(n_pnn=30, balanced=True, b_sight=3000, b_maxl=1500)

    links = oracle.get_links(cluster_name_for_GRN_unit=cluster_column, alpha=10, verbose_level=0)
    # p=0.001: significance threshold for TF-target connections
    links.filter_links(p=0.001, weight='coef_abs', threshold_number=2000)

    print(f'GRN constructed for {len(links.links_dict)} cell types')
    return oracle, links


def simulate_knockout(oracle, tf_name, n_propagation=3, n_neighbors=200, sigma_corr=0.05):
    '''Simulate TF knockout and compute cell state shifts.'''
    oracle.simulate_shift(perturb_condition={tf_name: 0.0}, n_propagation=n_propagation)
    oracle.estimate_transition_prob(n_neighbors=n_neighbors, knn_random=True, sampled_fraction=1)
    oracle.calculate_embedding_shift(sigma_corr=sigma_corr)

    shift = np.sqrt((oracle.adata.obsm['delta_embedding'] ** 2).sum(axis=1))
    oracle.adata.obs[f'{tf_name}_KO_shift'] = shift

    print(f'{tf_name} KO: mean shift = {shift.mean():.4f}, max = {shift.max():.4f}')
    return shift


def plot_perturbation(oracle, tf_name, save_prefix=None):
    '''Quiver plot and gradient plot for perturbation results.'''
    fig, axes = plt.subplots(1, 3, figsize=(24, 8))

    # Cell type reference
    sc.pl.embedding(oracle.adata, basis='umap', color='cell_type',
                    ax=axes[0], show=False, title='Cell types')

    # Quiver plot
    oracle.plot_quiver(ax=axes[1], scale=30, plot_whole_cells=True)
    axes[1].set_title(f'{tf_name} KO - predicted shifts')

    # Shift magnitude
    sc.pl.embedding(oracle.adata, basis='umap', color=f'{tf_name}_KO_shift',
                    cmap='Reds', ax=axes[2], show=False, title='Shift magnitude')

    plt.tight_layout()
    save_path = f'{save_prefix or tf_name}_ko_results.pdf'
    plt.savefig(save_path, bbox_inches='tight')
    plt.close()
    print(f'Saved {save_path}')


def screen_tfs(oracle, tf_list, n_propagation=3):
    '''Screen multiple TFs to find drivers of cell fate.'''
    results = {}
    for tf in tf_list:
        shift = simulate_knockout(oracle, tf, n_propagation=n_propagation)
        results[tf] = {
            'mean_shift': shift.mean(),
            'max_shift': shift.max(),
            'affected_cells_90pct': (shift > shift.quantile(0.9)).sum()
        }

    screen_df = pd.DataFrame(results).T.sort_values('mean_shift', ascending=False)
    screen_df.to_csv('tf_screen_results.csv')
    print('\nTF screen results:')
    print(screen_df)
    return screen_df


if __name__ == '__main__':
    # base_grn = build_base_grn('atac_peaks.bed')

    # adata = sc.read_h5ad('clustered.h5ad')
    # base_grn = pd.read_parquet('base_grn.parquet')
    # oracle, links = construct_grn(adata, base_grn)

    # Single TF knockout
    # shift = simulate_knockout(oracle, 'GATA1')
    # plot_perturbation(oracle, 'GATA1')

    # Screen multiple TFs
    # tfs = ['GATA1', 'SPI1', 'CEBPA', 'PAX5', 'TCF7', 'RUNX1']
    # screen_df = screen_tfs(oracle, tfs)

    print('Uncomment sections above to run with actual data')
