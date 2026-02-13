'''
Differential co-expression network analysis in Python.
Compares correlation networks between two conditions using Fisher z-transform.
'''
# Reference: matplotlib 3.8+, numpy 1.26+, pandas 2.2+, scipy 1.12+, statsmodels 0.14+ | Verify API if version differs

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def fisher_z_test(r1, n1, r2, n2):
    '''Fisher z-transform test for difference between two correlations.'''
    z1 = np.arctanh(np.clip(r1, -0.9999, 0.9999))
    z2 = np.arctanh(np.clip(r2, -0.9999, 0.9999))
    se = np.sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
    z_stat = (z1 - z2) / se
    p_value = 2 * stats.norm.sf(abs(z_stat))
    return z_stat, p_value


def differential_network_analysis(expr_cond1, expr_cond2, fdr_threshold=0.05, cor_threshold=0.3):
    '''
    Compare correlation networks between two conditions.

    Args:
        expr_cond1: DataFrame (samples x genes) for condition 1
        expr_cond2: DataFrame (samples x genes) for condition 2
        fdr_threshold: FDR cutoff for significance
        cor_threshold: minimum absolute correlation to classify as an edge
    '''
    genes = expr_cond1.columns.tolist()
    n1, n2 = len(expr_cond1), len(expr_cond2)
    cor1 = expr_cond1.corr()
    cor2 = expr_cond2.corr()

    results = []
    for i in range(len(genes)):
        for j in range(i + 1, len(genes)):
            r1, r2 = cor1.iloc[i, j], cor2.iloc[i, j]
            z_stat, p_val = fisher_z_test(r1, n1, r2, n2)
            results.append({
                'gene1': genes[i], 'gene2': genes[j],
                'cor_cond1': r1, 'cor_cond2': r2,
                'delta_cor': r2 - r1, 'z_stat': z_stat, 'p_value': p_val
            })

    df = pd.DataFrame(results)
    _, df['p_adj'], _, _ = multipletests(df['p_value'], method='fdr_bh')

    def classify(row):
        r1, r2 = abs(row['cor_cond1']), abs(row['cor_cond2'])
        s1, s2 = np.sign(row['cor_cond1']), np.sign(row['cor_cond2'])
        if r1 < cor_threshold and r2 >= cor_threshold:
            return 'gained'
        if r1 >= cor_threshold and r2 < cor_threshold:
            return 'lost'
        if r1 >= cor_threshold and r2 >= cor_threshold and s1 != s2:
            return 'reversed'
        return 'unchanged'

    df['edge_type'] = df.apply(classify, axis=1)
    return df


def find_rewired_hubs(rewired_df, top_n=20):
    '''Find genes involved in the most rewired connections.'''
    gene_counts = pd.concat([rewired_df['gene1'], rewired_df['gene2']]).value_counts()
    return gene_counts.head(top_n)


def plot_differential_network(rewired_df, top_genes=50, save='differential_network.pdf'):
    '''Visualize differential network with edge coloring by type.'''
    gene_counts = pd.concat([rewired_df['gene1'], rewired_df['gene2']]).value_counts()
    top = gene_counts.head(top_genes).index.tolist()
    subset = rewired_df[rewired_df['gene1'].isin(top) & rewired_df['gene2'].isin(top)]

    color_map = {'gained': '#2ca02c', 'lost': '#d62728', 'reversed': '#9467bd'}
    G = nx.Graph()

    for _, row in subset.iterrows():
        G.add_edge(row['gene1'], row['gene2'],
                    color=color_map[row['edge_type']], weight=abs(row['z_stat']))

    edge_colors = [G[u][v]['color'] for u, v in G.edges()]

    fig, ax = plt.subplots(figsize=(14, 14))
    pos = nx.spring_layout(G, k=2, seed=42)
    nx.draw_networkx(G, pos, edge_color=edge_colors, node_size=40,
                     font_size=6, width=0.5, alpha=0.7, ax=ax)

    legend_handles = [mpatches.Patch(color=c, label=l) for l, c in color_map.items()]
    ax.legend(handles=legend_handles, loc='upper left', fontsize=10)
    ax.set_title('Differential co-expression network', fontsize=14)

    plt.savefig(save, bbox_inches='tight')
    plt.close()
    print(f'Saved {save}')


if __name__ == '__main__':
    # Load expression data
    # expr_all = pd.read_csv('normalized_counts.csv', index_col=0).T
    # sample_info = pd.read_csv('sample_info.csv', index_col=0)

    # Split by condition
    # expr_ctrl = expr_all.loc[sample_info[sample_info['condition'] == 'control'].index]
    # expr_disease = expr_all.loc[sample_info[sample_info['condition'] == 'disease'].index]

    # Filter to top variable genes
    # gene_vars = expr_all.var()
    # top_genes = gene_vars.nlargest(3000).index.tolist()
    # expr_ctrl, expr_disease = expr_ctrl[top_genes], expr_disease[top_genes]

    # Run differential network analysis
    # diff_results = differential_network_analysis(expr_ctrl, expr_disease)
    # significant = diff_results[diff_results['p_adj'] < 0.05]
    # rewired = significant[significant['edge_type'] != 'unchanged']

    # print(f'Significant pairs: {len(significant)}')
    # print(significant['edge_type'].value_counts())

    # Hub genes
    # hubs = find_rewired_hubs(rewired)
    # print(f'\nTop rewired hub genes:\n{hubs}')

    # Visualize
    # plot_differential_network(rewired)
    # rewired.to_csv('rewired_edges.csv', index=False)

    print('Uncomment sections above to run with actual data')
