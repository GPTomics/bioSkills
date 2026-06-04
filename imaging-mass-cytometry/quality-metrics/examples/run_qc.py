'''IMC channel QC: count-aware cell-level SNR, empty-channel test, sample-of-origin check.

Replaces the fluorescence-style "signal mean / background std" SNR with the
decision-relevant metric: a two-component Gaussian mixture on per-cell COUNTS (a marker
fails when positive and negative do not separate, regardless of brightness), plus the
operational "is this antibody dead" test (statistically indistinguishable from an empty
channel), plus the dominant batch failure (cells clustering by sample-of-origin).
'''
# Reference: numpy 1.26+, scikit-learn 1.4+, pandas 2.2+ | Verify API if version differs
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture

rng = np.random.default_rng(0)
n_cells = 4000

# simulate per-cell ion COUNTS (Poisson): a real marker (two populations), a dead antibody
# (one population near the floor), and an empty/background channel
real = np.where(rng.random(n_cells) < 0.3, rng.poisson(12, n_cells), rng.poisson(1, n_cells))
dead = rng.poisson(1, n_cells)
empty = rng.poisson(1, n_cells)

def cell_snr(counts):
    # two-component mixture on raw counts, reported as context: on toy Poisson data a GMM can
    # split pure noise into a spurious "positive" component, which is exactly why the
    # empty-channel comparison below -- not the raw SNR ratio -- is the decisive test
    gm = GaussianMixture(n_components=2, random_state=0).fit(counts.reshape(-1, 1).astype(float))
    hi, lo = np.sort(gm.means_.ravel())[::-1]
    return hi / lo if lo > 0 else np.inf

def signal_above_background(counts, empty_counts, q=95, tol=2.0):
    # operational "did this antibody work" test: is the POSITIVE tail above the empty-channel
    # floor? A dim-but-real marker carries its signal in the tail; a failed channel's tail sits
    # at the empty floor. Comparing medians would wrongly fail a sparse-but-real marker whose
    # median is the negative floor -- exactly the dim != failed trap the skill warns about.
    return np.percentile(counts, q) - np.percentile(empty_counts, q) > tol

for name, counts in [('real_marker', real), ('dead_antibody', dead)]:
    snr = cell_snr(counts)
    works = signal_above_background(counts, empty)
    verdict = 'KEEP' if works else 'DROP (failed: signal ~ empty channel)'
    print(f'{name}: cell-SNR={snr:.2f}  signal_above_background={works}  -> {verdict}')

# dominant batch failure: cells separating by which slide they came from, not phenotype.
# Here a fabricated per-slide offset stands in for FFPE/ischemia/lot variation; in a real
# pipeline, UMAP the cell x marker matrix and color by patient/slide/day/antibody_lot.
slide = rng.integers(0, 3, n_cells)
marker_with_batch = real + slide * 5.0
by_slide = pd.DataFrame({'slide': slide, 'value': marker_with_batch}).groupby('slide')['value'].median()
print('\nPer-slide median (large spread = sample-of-origin batch effect, gate before analysis):')
print(by_slide.round(1).to_string())
