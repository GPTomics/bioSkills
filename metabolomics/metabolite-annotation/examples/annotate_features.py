'''Spectral library matching with a confidence-level call on toy MS/MS spectra.

Demonstrates the skill's central claim: a match score ranks candidates, it does
not identify a compound. The level returned is a hypothesis, never an identification.
Self-contained -- builds synthetic Spectrum objects, no external database needed.
'''
# Reference: matchms 0.33+ | Verify API if version differs
import numpy as np
from matchms import Spectrum, calculate_scores
from matchms.filtering import normalize_intensities, add_precursor_mz

try:
    from matchms.similarity import ModifiedCosineGreedy as ModifiedCosine  # matchms 0.33+
except ImportError:
    from matchms.similarity import ModifiedCosine  # matchms <= 0.32

SCORE_FLOOR = 0.7   # GNPS default; below this, alignment is dominated by noise
MATCH_FLOOR = 6     # matched-peak floor; a high score on few peaks is meaningless


def make_spectrum(name, precursor, mz, intensities):
    spectrum = Spectrum(mz=np.array(mz, dtype=float), intensities=np.array(intensities, dtype=float),
                        metadata={'compound_name': name, 'precursor_mz': precursor})
    spectrum = add_precursor_mz(spectrum)  # ModifiedCosine returns zeros without this
    return normalize_intensities(spectrum)


references = [
    make_spectrum('hippuric_acid', 180.0655, [65.04, 77.04, 91.05, 105.03, 134.06, 162.05, 180.06],
                  [0.2, 0.4, 0.25, 1.0, 0.6, 0.3, 0.15]),
    make_spectrum('phenylacetylglycine', 194.0812, [65.04, 76.03, 91.05, 118.06, 148.08, 176.07, 194.08],
                  [0.15, 0.3, 1.0, 0.5, 0.4, 0.2, 0.18]),
]

# Two queries: one that strongly matches hippuric acid (many shared peaks), one
# low-information spectrum (two generic peaks) that scores high on too few peaks --
# the classic over-claim trap the matched-peak floor is designed to catch.
queries = [
    make_spectrum('query_strong', 180.0655, [65.04, 77.04, 91.05, 105.03, 134.06, 162.05, 180.06],
                  [0.21, 0.38, 0.24, 1.0, 0.58, 0.31, 0.16]),
    make_spectrum('query_promiscuous', 180.0655, [91.05, 105.03], [0.9, 1.0]),
]

scores = calculate_scores(references, queries, ModifiedCosine(tolerance=0.01))


def assign_level(score, matches, has_standard=False):
    if has_standard:
        return 1
    if score >= SCORE_FLOOR and matches >= MATCH_FLOOR:
        return '2a'  # reference library match, no in-house standard
    if score >= SCORE_FLOOR:
        return 3     # high score but too few peaks -> candidate only, isomers unresolved
    return 5         # no defensible structural evidence


def best_hit(query):
    pairs = scores.scores_by_query(query)
    score_field, match_field = pairs[0][1].dtype.names  # names are class-prefixed and version-dependent
    ref, hit = max(pairs, key=lambda pair: pair[1][score_field])
    return ref, float(hit[score_field]), int(hit[match_field])


for query in queries:
    ref, score, matches = best_hit(query)
    level = assign_level(score, matches)
    print(f"{query.get('compound_name'):>18} -> {ref.get('compound_name'):<20} score={score:.2f} matches={matches} level={level}")
