"""CLR-transform a species-by-sample abundance matrix the compositionally-correct way.

A shotgun abundance table is a composition: counts sum to a sequencer-fixed total, so raw
proportions carry only relative information and correlations of proportions are biased. Replace
zeros (log of zero is undefined), CLR-transform to remove closure, then use Aitchison distance
(Euclidean on CLR) for ordination - never Pearson or Bray-Curtis on raw proportions.
"""
# Reference: pandas 2.2+, scikit-bio 0.6+, numpy 1.26+ | Verify API if version differs
import numpy as np
import pandas as pd
from skbio.stats.composition import clr, multi_replace
from scipy.spatial.distance import pdist, squareform

# scikit-bio renamed multiplicative_replacement -> multi_replace in 0.6.0; introspect if it differs.


def clr_transform(counts):
    """counts: taxa x samples DataFrame of read counts (e.g. Bracken new_est_reads)."""
    mat = counts.T.values.astype(float)
    proportions = mat / mat.sum(axis=1, keepdims=True)
    no_zero = multi_replace(proportions)
    clr_mat = clr(no_zero)
    return pd.DataFrame(clr_mat, index=counts.columns, columns=counts.index)


def aitchison_distance(clr_df):
    """Euclidean distance on CLR coordinates = Aitchison distance on the original composition."""
    dist = squareform(pdist(clr_df.values, metric='euclidean'))
    return pd.DataFrame(dist, index=clr_df.index, columns=clr_df.index)


def redistribution_flag(bracken_df, added_to_assigned_ratio=5.0):
    """Flag species whose abundance is mostly redistributed - possible artifacts of DB-absent taxa.

    added_to_assigned_ratio default 5.0 is a heuristic here, not a Bracken parameter: a species with
    >5x more added than directly-assigned (kraken_assigned_reads) reads is suspicious (Bracken
    redistributed into a database-present relative of an absent organism).
    """
    assigned = bracken_df['kraken_assigned_reads'].replace(0, 1)
    ratio = bracken_df['added_reads'] / assigned
    return bracken_df.loc[ratio > added_to_assigned_ratio, 'name'].tolist()


if __name__ == '__main__':
    rng = np.random.default_rng(0)
    taxa = [f'species_{i}' for i in range(8)]
    samples = [f'sample_{j}' for j in range(6)]
    counts = pd.DataFrame(rng.poisson(lam=rng.uniform(1, 500, size=(8, 6))), index=taxa, columns=samples)
    counts.iloc[5:, :3] = 0  # structural-looking zeros to exercise replacement

    clr_df = clr_transform(counts)
    dist = aitchison_distance(clr_df)

    print('CLR coordinates (samples x taxa):')
    print(clr_df.round(2))
    print('\nAitchison distance between samples:')
    print(dist.round(2))
    assert np.allclose(clr_df.values.sum(axis=1), 0, atol=1e-9), 'CLR rows must sum to zero'
    print('\nCLR rows sum to zero (closure removed): OK')
