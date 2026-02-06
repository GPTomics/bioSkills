'''
Visualize pySCENIC results: regulon heatmaps, RSS plots, and UMAP overlays.

Requires:
- clustered.h5ad: preprocessed scRNA-seq with UMAP and cell type labels
- auc_matrix.csv: AUCell scores from pySCENIC
- rss_scores.csv: regulon specificity scores
'''

import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pyscenic.binarization import binarize


def plot_regulon_heatmap(auc_mtx, cell_types, top_n=5, save='regulon_heatmap.pdf'):
    '''Heatmap of mean regulon activity per cell type, showing top N per type.'''
    mean_activity = auc_mtx.groupby(cell_types).mean()

    # Select top regulons per cell type by mean activity
    top_regulons = set()
    for ct in mean_activity.index:
        top = mean_activity.loc[ct].nlargest(top_n).index.tolist()
        top_regulons.update(top)
    top_regulons = sorted(top_regulons)

    g = sns.clustermap(mean_activity[top_regulons].T, cmap='viridis',
                       figsize=(12, max(6, len(top_regulons) * 0.3)),
                       z_score=0, xticklabels=True, yticklabels=True)
    g.fig.suptitle('Regulon activity by cell type (z-scored)', y=1.02)
    plt.savefig(save, bbox_inches='tight')
    plt.close()
    print(f'Saved {save}')


def plot_rss_panel(rss, cell_types_to_show=None, top_n=5, save='rss_panel.pdf'):
    '''Panel of RSS rank plots per cell type.'''
    if cell_types_to_show is None:
        cell_types_to_show = rss.columns.tolist()

    n_types = len(cell_types_to_show)
    ncols = min(4, n_types)
    nrows = (n_types + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows))
    axes = np.array(axes).flatten()

    for i, ct in enumerate(cell_types_to_show):
        ax = axes[i]
        scores = rss[ct].sort_values(ascending=False)
        ax.plot(range(len(scores)), scores.values, '.', color='steelblue', markersize=3)

        for j, (reg, score) in enumerate(scores.head(top_n).items()):
            ax.annotate(reg, (j, score), fontsize=7, ha='left', va='bottom')

        ax.set_title(ct, fontsize=10)
        ax.set_xlabel('Regulon rank')
        ax.set_ylabel('RSS')

    for i in range(n_types, len(axes)):
        axes[i].set_visible(False)

    plt.tight_layout()
    plt.savefig(save, bbox_inches='tight')
    plt.close()
    print(f'Saved {save}')


def plot_regulon_umap(adata, auc_mtx, regulons_to_plot, save='regulon_umap.pdf'):
    '''Overlay regulon activity on UMAP embedding.'''
    for reg in regulons_to_plot:
        if reg in auc_mtx.columns:
            adata.obs[reg] = auc_mtx.loc[adata.obs_names, reg].values

    fig = sc.pl.umap(adata, color=regulons_to_plot, cmap='viridis',
                     ncols=3, return_fig=True)
    fig.savefig(save, bbox_inches='tight')
    plt.close()
    print(f'Saved {save}')


def plot_binary_activity(auc_mtx, cell_types, save='binary_activity.pdf'):
    '''Heatmap of fraction of cells with active regulon per cell type.'''
    binary_mtx, thresholds = binarize(auc_mtx)
    fraction_active = binary_mtx.groupby(cell_types).mean()

    # Keep regulons active in at least one cell type (>10% of cells)
    active_regs = fraction_active.columns[fraction_active.max() > 0.1]
    fraction_active = fraction_active[active_regs]

    g = sns.clustermap(fraction_active.T, cmap='YlOrRd',
                       figsize=(12, max(6, len(active_regs) * 0.25)),
                       vmin=0, vmax=1, xticklabels=True, yticklabels=True)
    g.fig.suptitle('Fraction of cells with active regulon', y=1.02)
    plt.savefig(save, bbox_inches='tight')
    plt.close()
    print(f'Saved {save}')


if __name__ == '__main__':
    # adata = sc.read_h5ad('clustered.h5ad')
    # auc_mtx = pd.read_csv('auc_matrix.csv', index_col=0)
    # rss = pd.read_csv('rss_scores.csv', index_col=0)
    # cell_types = adata.obs['cell_type']

    # plot_regulon_heatmap(auc_mtx, cell_types)
    # plot_rss_panel(rss)
    # plot_regulon_umap(adata, auc_mtx, ['CEBPB(+)', 'SPI1(+)', 'PAX5(+)'])
    # plot_binary_activity(auc_mtx, cell_types)

    print('Uncomment sections above to run with actual pySCENIC results')
