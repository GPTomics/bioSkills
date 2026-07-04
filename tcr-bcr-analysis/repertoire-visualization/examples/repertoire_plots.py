'''Repertoire visualization recipes for TCR/BCR data.

Each function draws one figure type and returns a matplotlib Figure. The recipes
encode the interpretation, not just the drawing: spectratype shape reads as
poly- vs oligoclonal; diversity is compared by rarefaction (not a raw bar);
overlap uses a depth-robust metric; the network structure is threshold-dependent.
Runs end-to-end on a synthetic clonotype DataFrame and leaves no output files.
'''
# Reference: matplotlib 3.8+ seaborn 0.13+ | Verify API if version differs

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

STRATA_EDGES = [0, 1e-4, 1e-3, 1e-2, 1e-1, 1.0]  # frequency bins; frequency-based makes clonal-space depth-robust
STRATA_LABELS = ['Rare', 'Small', 'Medium', 'Large', 'Hyperexpanded']


def plot_spectratype(clone_df, length_col='cdr3_length'):
    '''CDR3-length histogram. Gaussian = polyclonal/naive; spikes = clonal expansion.

    Weight by frequency to surface expansions; by unique clonotypes to show diversity.
    '''
    fig, ax = plt.subplots(figsize=(9, 5))
    bins = range(int(clone_df[length_col].min()), int(clone_df[length_col].max()) + 2)
    ax.hist(clone_df[length_col], bins=bins, weights=clone_df['frequency'], color='steelblue')
    ax.set_xlabel('CDR3 length (aa)')
    ax.set_ylabel('Frequency')
    ax.set_title('CDR3 spectratype (frequency-weighted)')
    fig.tight_layout()
    return fig


def plot_clonal_space(clone_df, sample_col='sample'):
    '''Stacked bar of clone-size strata (clonal-space homeostasis).'''
    clone_df = clone_df.copy()
    clone_df['stratum'] = pd.cut(clone_df['frequency'], bins=STRATA_EDGES, labels=STRATA_LABELS)
    space = clone_df.groupby([sample_col, 'stratum'], observed=True)['frequency'].sum().unstack(fill_value=0)
    space = space.reindex(columns=STRATA_LABELS, fill_value=0)
    fig, ax = plt.subplots(figsize=(8, 5))
    space.plot(kind='bar', stacked=True, ax=ax, colormap='viridis')
    ax.set_ylabel('Fraction of repertoire')
    ax.set_title('Clonal-space homeostasis')
    ax.legend(title='Stratum', bbox_to_anchor=(1.02, 1), loc='upper left')
    fig.tight_layout()
    return fig


def plot_clone_tracking(clone_df, top_n=10, clone_col='cdr3_aa', time_col='timepoint'):
    '''Line plot of the top clones across timepoints.

    Downsample timepoints to common depth before calling contraction; an absent
    clone is often a sampling zero rather than a true loss.
    '''
    top = clone_df.groupby(clone_col)['frequency'].sum().nlargest(top_n).index
    fig, ax = plt.subplots(figsize=(9, 5))
    for clone in top:
        d = clone_df[clone_df[clone_col] == clone].sort_values(time_col)
        ax.plot(d[time_col], d['frequency'], marker='o', label=clone[:12])
    ax.set_xlabel('Timepoint')
    ax.set_ylabel('Clone frequency')
    ax.set_title(f'Top {top_n} clone dynamics')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=7)
    fig.tight_layout()
    return fig


def rarefaction_curve(counts, depths, reps=20, rng=None):
    '''Interpolated observed richness vs sampling depth by multinomial resampling.

    The correct way to compare diversity across unequal depth: read curves at a
    shared x. Interpolation only here - do not extrapolate past observed depth.
    '''
    rng = rng or np.random.default_rng(0)
    counts = np.asarray(counts, dtype=float)
    p = counts / counts.sum()
    total = int(counts.sum())
    richness = []
    for m in depths:
        if m > total:
            richness.append(np.nan)
            continue
        obs = [np.count_nonzero(rng.multinomial(int(m), p)) for _ in range(reps)]
        richness.append(float(np.mean(obs)))
    return richness


def plot_rarefaction(sample_counts, rng=None):
    '''Rarefaction curves per sample, read at a common x rather than a raw bar.'''
    fig, ax = plt.subplots(figsize=(8, 5))
    min_depth = min(int(np.sum(c)) for c in sample_counts.values())
    depths = np.linspace(10, min_depth, 12).astype(int)  # common x-range = shallowest sample
    for name, counts in sample_counts.items():
        ax.plot(depths, rarefaction_curve(counts, depths, rng=rng), marker='.', label=name)
    ax.axvline(min_depth, ls='--', color='grey', lw=1, label='common depth')
    ax.set_xlabel('Sampled reads')
    ax.set_ylabel('Observed clonotypes')
    ax.set_title('Rarefaction (compare diversity here, not with a raw bar)')
    ax.legend(fontsize=8)
    fig.tight_layout()
    return fig


def plot_overlap_heatmap(overlap_matrix, metric='Morisita-Horn'):
    '''Pairwise overlap heatmap. Morisita-Horn is depth-robust; Jaccard is not.'''
    fig, ax = plt.subplots(figsize=(7, 6))
    sns.heatmap(overlap_matrix, annot=True, fmt='.2f', cmap='YlOrRd', vmin=0, vmax=1, square=True, ax=ax)
    ax.set_title(f'Repertoire overlap ({metric})')
    fig.tight_layout()
    return fig


def plot_vj_heatmap(clone_df):
    '''V-by-J pairing frequency heatmap (quantitative alternative to a chord).'''
    vj = clone_df.pivot_table(index='v_gene', columns='j_gene', values='frequency', aggfunc='sum', fill_value=0)
    fig, ax = plt.subplots(figsize=(7, 5))
    sns.heatmap(vj, cmap='viridis', ax=ax)
    ax.set_title('V-J pairing frequency')
    fig.tight_layout()
    return fig


def _levenshtein(a, b):
    if a == b:
        return 0
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def plot_similarity_network(clone_df, max_norm_dist=0.2, clone_col='cdr3_aa'):
    '''Clonotype-similarity network. Structure depends entirely on the threshold.'''
    clones = list(clone_df[clone_col].unique())
    n = len(clones)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    pos = {c: (np.cos(a), np.sin(a)) for c, a in zip(clones, angles)}
    fig, ax = plt.subplots(figsize=(7, 7))
    for i, a in enumerate(clones):
        for b in clones[i + 1:]:
            d = _levenshtein(a, b) / max(len(a), len(b))
            if d <= max_norm_dist:
                (x0, y0), (x1, y1) = pos[a], pos[b]
                ax.plot([x0, x1], [y0, y1], color='grey', lw=0.8, alpha=0.6)
    xs, ys = zip(*pos.values())
    ax.scatter(xs, ys, s=60, color='steelblue', zorder=3)
    ax.set_title(f'CDR3 similarity network (norm dist <= {max_norm_dist})')
    ax.axis('off')
    fig.tight_layout()
    return fig


def _synthetic_repertoire(rng):
    aa = list('ACDEFGHIKLMNPQRSTVWY')
    v_genes = ['TRBV5-1', 'TRBV6-1', 'TRBV7-2', 'TRBV20-1']
    j_genes = ['TRBJ1-1', 'TRBJ2-1', 'TRBJ2-7']
    rows = []
    for sample, tp in [('S1', 1), ('S2', 2), ('S3', 3)]:
        n = 40
        counts = (rng.zipf(1.8, n)).astype(float)  # long-tailed abundance with a few expanded clones
        counts = np.clip(counts, 1, 500)
        freq = counts / counts.sum()
        lengths = np.clip(rng.normal(15, 2, n).round().astype(int), 10, 22)
        for k in range(n):
            cdr3 = 'C' + ''.join(rng.choice(aa, size=lengths[k] - 2)) + 'F'
            rows.append({'sample': sample, 'timepoint': tp, 'cdr3_aa': cdr3,
                         'cdr3_length': lengths[k], 'frequency': freq[k], 'count': counts[k],
                         'v_gene': rng.choice(v_genes), 'j_gene': rng.choice(j_genes)})
    return pd.DataFrame(rows)


if __name__ == '__main__':
    import os
    import shutil
    import tempfile

    rng = np.random.default_rng(0)
    df = _synthetic_repertoire(rng)

    overlap = pd.DataFrame([[1.0, 0.32, 0.11], [0.32, 1.0, 0.28], [0.11, 0.28, 1.0]],
                           index=['S1', 'S2', 'S3'], columns=['S1', 'S2', 'S3'])
    sample_counts = {s: df[df['sample'] == s]['count'].values for s in ['S1', 'S2', 'S3']}

    figs = {
        'spectratype': plot_spectratype(df),
        'clonal_space': plot_clonal_space(df),
        'clone_tracking': plot_clone_tracking(df),
        'rarefaction': plot_rarefaction(sample_counts, rng=rng),
        'overlap': plot_overlap_heatmap(overlap),
        'vj_heatmap': plot_vj_heatmap(df),
        'network': plot_similarity_network(df),
    }

    tmp = tempfile.mkdtemp(prefix='repviz_')
    for name, fig in figs.items():
        fig.savefig(os.path.join(tmp, f'{name}.png'), dpi=90)
        plt.close(fig)
    shutil.rmtree(tmp)  # leave no output files in the repo
    print(f'Rendered {len(figs)} figures on {len(df)} clonotypes; temp output cleaned.')
