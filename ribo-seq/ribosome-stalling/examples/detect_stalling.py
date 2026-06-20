'''Detect ribosome pause sites with a local-relative pause score.

Pause score = codon occupancy / gene-mean occupancy, computed per transcript.
This local-relative metric (not a global z-score across genes) is the field
standard. Valid only on flash-frozen, no-drug data: cycloheximide pre-treatment
redistributes density and fabricates pauses (Hussmann 2015). Occupancy is
assigned to the A-site (A-site offset = P-site + 3) for decoding/tRNA effects.
Real per-CDS occupancy comes from plastid cds.get_counts(BAMGenomeArray); this
script demonstrates the scoring on a simulated occupancy vector.
'''
# Reference: numpy 1.26+, biopython 1.83+ | Verify API if version differs

import numpy as np

MIN_TOTAL = 500       # per-gene footprint floor; below this, pause scores are noise
SCORE_THRESHOLD = 5.0  # fold-over-gene-mean to call a pause (tune per dataset)


def pause_scores(per_codon_occupancy, min_total=MIN_TOTAL, threshold=SCORE_THRESHOLD):
    '''Codon pause scores = occupancy / gene-mean occupancy, above a coverage floor.'''
    pauses = []
    for tx, occ in per_codon_occupancy.items():
        occ = np.asarray(occ, dtype=float)
        if occ.sum() < min_total:
            continue
        mean = occ.mean()
        if mean == 0:
            continue
        scores = occ / mean
        for pos in np.where(scores > threshold)[0]:
            pauses.append({'transcript': tx, 'codon': int(pos),
                           'pause_score': float(scores[pos])})
    return pauses


def codon_occupancy(per_codon_occupancy, codon_seqs):
    '''Per-codon-type occupancy, normalizing each gene to its own mean BEFORE pooling.

    Mean-of-ratios (not ratio-of-means) so highly expressed genes do not dominate.
    codon_seqs maps transcript -> list of codon strings aligned to the occupancy vector.
    '''
    by_codon = {}
    for tx, occ in per_codon_occupancy.items():
        occ = np.asarray(occ, dtype=float)
        mean = occ.mean()
        if mean == 0:
            continue
        norm = occ / mean
        for codon, value in zip(codon_seqs.get(tx, []), norm):
            by_codon.setdefault(codon, []).append(value)
    return {c: float(np.mean(v)) for c, v in by_codon.items() if len(v) >= 100}


def _simulate(n_tx=50, n_codons=200, depth=8, seed=0):
    '''Simulate per-codon A-site occupancy with a planted pause in each transcript.'''
    rng = np.random.default_rng(seed)
    occ, seqs = {}, {}
    codon_pool = ['CCG', 'AAA', 'GAA', 'CTG', 'GCC', 'ATG']
    for t in range(n_tx):
        v = rng.poisson(depth, n_codons).astype(float)
        v[rng.integers(20, n_codons - 20)] += depth * 12   # planted pause
        occ[f'tx{t}'] = v
        seqs[f'tx{t}'] = list(rng.choice(codon_pool, n_codons))
    return occ, seqs


if __name__ == '__main__':
    print('Ribosome stalling: local pause score (simulated A-site occupancy)')
    occ, seqs = _simulate()
    pauses = pause_scores(occ)
    print(f'  pauses called: {len(pauses)} (expect ~1 per transcript)')
    if pauses:
        top = max(pauses, key=lambda p: p['pause_score'])
        print(f"  strongest: {top['transcript']} codon {top['codon']} "
              f"score={top['pause_score']:.1f}")
    means = codon_occupancy(occ, seqs)
    print(f'  per-codon-type occupancy estimated for {len(means)} codons')
