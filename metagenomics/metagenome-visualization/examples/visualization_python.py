"""Compositionally coherent visualization of a shotgun profiler table.

The ordination here is CLR-PCA (Aitchison-PCA), NOT StandardScaler-then-PCA on raw relative abundance,
which is compositionally incoherent. The stacked bar labels how many taxa are hidden in "Other" - a
relative-abundance bar shows the relative race, never the absolute community.
"""
# Reference: pandas 2.2+, scikit-bio 0.6+, scikit-learn 1.4+, matplotlib 3.8+ | Verify API if version differs
import os
import sys
import tempfile
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from skbio.stats.composition import clr, multi_replace   # renamed from multiplicative_replacement in skbio 0.6
from sklearn.decomposition import PCA


def load_species(path):
    ab = pd.read_csv(path, sep='\t', index_col=0)
    species = ab[ab.index.str.contains(r'\|s__') & ~ab.index.str.contains(r'\|t__')].copy()
    species.index = species.index.str.split('|').str[-1].str.replace('s__', '')
    return species


def clr_pca(species, n_components=2):
    proportions = species.T.values / species.T.values.sum(axis=1, keepdims=True)
    clr_mat = clr(multi_replace(proportions))
    pca = PCA(n_components=n_components).fit(clr_mat)
    return pca, pca.transform(clr_mat)


def honest_stacked_bar(species, top_n, out_path):
    top = species.sum(axis=1).nlargest(top_n).index
    plotted = species.loc[top].copy()
    hidden = species.drop(top)
    plotted.loc['Other'] = hidden.sum()
    hidden_pct = 100 * hidden.sum().sum() / species.sum().sum()
    fig, ax = plt.subplots(figsize=(10, 6))
    plotted.T.plot(kind='bar', stacked=True, ax=ax, colormap='tab20')
    ax.set_ylabel('Relative abundance (%)')
    ax.set_title(f'Composition (Other hides {len(hidden)} taxa, {hidden_pct:.1f}% of signal; relative, not absolute)')
    plt.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return len(hidden), hidden_pct


if __name__ == '__main__':
    path = sys.argv[1] if len(sys.argv) > 1 else None
    if path:
        species = load_species(path)
    else:
        rng = np.random.default_rng(0)
        taxa = [f'species_{i}' for i in range(25)]
        samples = [f's{j}' for j in range(8)]
        species = pd.DataFrame(rng.poisson(rng.uniform(1, 300, size=(25, 8))), index=taxa, columns=samples).astype(float)

    out_path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(tempfile.gettempdir(), 'stacked_bar.png')
    pca, coords = clr_pca(species)
    n_hidden, pct = honest_stacked_bar(species, 10, out_path)
    print(f'CLR-PCA variance explained: {(pca.explained_variance_ratio_[:2] * 100).round(1)}%')
    print(f'Stacked bar Other hides {n_hidden} taxa ({pct:.1f}% of signal)')
    print(f'Wrote {out_path} (CLR-PCA coords computed; loadings interpretable as taxa)')
