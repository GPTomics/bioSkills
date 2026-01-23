'''Repertoire visualization examples'''

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_diversity_comparison(diversity_df, metric='shannon', groupby='condition'):
    '''Compare diversity metrics between groups

    Args:
        diversity_df: DataFrame with diversity metrics per sample
        metric: Diversity metric to plot (shannon, gini, chao1, d50)
        groupby: Column to group samples by
    '''
    fig, ax = plt.subplots(figsize=(8, 6))

    sns.boxplot(data=diversity_df, x=groupby, y=metric, ax=ax, palette='Set2')
    sns.stripplot(data=diversity_df, x=groupby, y=metric, ax=ax,
                  color='black', alpha=0.5, size=8)

    ax.set_ylabel(f'{metric.replace("_", " ").title()}')
    ax.set_xlabel('')
    plt.tight_layout()
    plt.savefig(f'diversity_{metric}.pdf', dpi=300)
    plt.close()

def plot_spectratype(clone_df, sample_col='sample', color_palette='tab10'):
    '''Plot CDR3 length distribution (spectratype)

    Normal repertoire shows Gaussian distribution centered ~45 nt (15 aa)
    Skewed or multi-modal distributions may indicate oligoclonality
    '''
    clone_df = clone_df.copy()
    clone_df['cdr3_length'] = clone_df['cdr3_aa'].str.len()

    fig, ax = plt.subplots(figsize=(10, 6))

    samples = clone_df[sample_col].unique()
    colors = plt.cm.get_cmap(color_palette, len(samples))

    for i, sample in enumerate(samples):
        data = clone_df[clone_df[sample_col] == sample]
        # Weight by frequency for accurate representation
        weights = data['frequency'].values
        ax.hist(data['cdr3_length'], bins=range(8, 25), weights=weights,
                alpha=0.5, label=sample, color=colors(i), density=True)

    ax.set_xlabel('CDR3 Length (amino acids)')
    ax.set_ylabel('Frequency')
    ax.set_title('CDR3 Length Distribution')
    ax.legend()
    ax.set_xlim(8, 24)
    plt.tight_layout()
    plt.savefig('spectratype.pdf', dpi=300)
    plt.close()

def plot_clone_tracking(clone_df, timepoint_col='timepoint', top_n=10):
    '''Track top clones across timepoints

    Useful for vaccination, infection response, or tumor monitoring studies
    '''
    # Sum frequencies across timepoints to find top clones
    total_freq = clone_df.groupby('cdr3_aa')['frequency'].sum()
    top_clones = total_freq.nlargest(top_n).index.tolist()

    fig, ax = plt.subplots(figsize=(10, 6))

    for clone in top_clones:
        clone_data = clone_df[clone_df['cdr3_aa'] == clone].sort_values(timepoint_col)
        # Truncate label for readability
        label = clone[:15] + '...' if len(clone) > 15 else clone
        ax.plot(clone_data[timepoint_col], clone_data['frequency'],
                marker='o', linewidth=2, markersize=8, label=label)

    ax.set_xlabel('Timepoint')
    ax.set_ylabel('Clone Frequency')
    ax.set_title(f'Top {top_n} Clone Dynamics')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    ax.set_ylim(bottom=0)
    plt.tight_layout()
    plt.savefig('clone_tracking.pdf', dpi=300, bbox_inches='tight')
    plt.close()

def plot_overlap_heatmap(overlap_matrix, title='Repertoire Overlap'):
    '''Plot pairwise repertoire overlap as heatmap

    Args:
        overlap_matrix: Square DataFrame with sample pairs
        title: Plot title
    '''
    fig, ax = plt.subplots(figsize=(10, 8))

    # Use diverging colormap centered at median
    sns.heatmap(
        overlap_matrix,
        annot=True,
        fmt='.2f',
        cmap='YlOrRd',
        ax=ax,
        vmin=0,
        vmax=1,
        square=True
    )

    ax.set_title(title)
    plt.tight_layout()
    plt.savefig('overlap_heatmap.pdf', dpi=300)
    plt.close()

def plot_vgene_usage(clone_df, gene_col='v_gene', sample_col='sample', top_n=15):
    '''Plot V gene usage comparison'''

    # Aggregate by V gene
    usage = clone_df.groupby([sample_col, gene_col])['frequency'].sum().unstack(fill_value=0)

    # Keep top N genes by total usage
    top_genes = usage.sum().nlargest(top_n).index
    usage = usage[top_genes]

    fig, ax = plt.subplots(figsize=(12, 6))

    usage.T.plot(kind='bar', ax=ax, width=0.8)

    ax.set_xlabel('V Gene')
    ax.set_ylabel('Total Frequency')
    ax.set_title('V Gene Usage by Sample')
    ax.legend(title='Sample', bbox_to_anchor=(1.05, 1))
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig('vgene_usage.pdf', dpi=300, bbox_inches='tight')
    plt.close()

# Example usage
if __name__ == '__main__':
    # Create example data
    clone_df = pd.DataFrame({
        'cdr3_aa': ['CASSF', 'CASSP', 'CASSK', 'CASSL'] * 3,
        'frequency': [0.1, 0.05, 0.03, 0.02] * 3,
        'sample': ['S1'] * 4 + ['S2'] * 4 + ['S3'] * 4,
        'timepoint': [1] * 4 + [2] * 4 + [3] * 4,
        'v_gene': ['TRBV5-1', 'TRBV6-1', 'TRBV7-2', 'TRBV5-1'] * 3
    })

    plot_spectratype(clone_df)
    plot_vgene_usage(clone_df)
    print('Plots saved!')
