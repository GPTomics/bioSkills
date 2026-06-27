'''Test the piRNA ping-pong signature and a tRF end-precision metric.

The ping-pong signature is a sharp excess of 10-nt 5'-5' overlaps between sense and
antisense reads (slicer-driven secondary piRNA biogenesis). End precision distinguishes
a processed small RNA (sharp 5' terminus) from random degradation (smeared ends).
'''
# Reference: numpy 1.26+, pandas 2.2+ | Verify API if version differs

import numpy as np
from collections import Counter


def ping_pong_zscore(plus_5p, minus_5p, max_overlap=30):
    '''Score 5'-5' overlaps between strands; a 10-nt spike is the ping-pong signature.

    plus_5p / minus_5p map a genomic 5' coordinate to a read count, per strand.
    A sense read at i and an antisense read whose 5' end is at i+overlap-1 share an
    'overlap'-nt 5' overlap. A z >> 0 at overlap 10 is evidence of an active pathway.
    '''
    hist = np.zeros(max_overlap + 1)
    for pos, n in plus_5p.items():
        for overlap in range(1, max_overlap + 1):
            partner = pos + overlap - 1
            if partner in minus_5p:
                hist[overlap] += n * minus_5p[partner]
    others = np.concatenate([hist[1:10], hist[11:]])
    z10 = (hist[10] - others.mean()) / (others.std() + 1e-9)
    return hist, z10


def end_precision(read_5p_positions):
    '''Fraction of reads sharing the modal 5' end (near 1.0 = precise, low = decay-like).'''
    c = Counter(read_5p_positions)
    return max(c.values()) / sum(c.values())


def _simulate_pingpong(n_loci=200, depth=20, seed=0):
    rng = np.random.default_rng(seed)
    plus, minus = {}, {}
    for _ in range(n_loci):
        i = int(rng.integers(0, 100000))
        plus[i] = plus.get(i, 0) + depth
        partner = i + 9          # 10-nt 5'-5' overlap -> antisense 5' at i + 10 - 1
        minus[partner] = minus.get(partner, 0) + depth
    return plus, minus


if __name__ == '__main__':
    plus, minus = _simulate_pingpong()
    hist, z10 = ping_pong_zscore(plus, minus)
    print(f'ping-pong z-score at 10 nt overlap: {z10:.1f} (>> 0 indicates an active pathway)')

    processed = [100] * 45 + [101] * 5
    decay = list(range(100, 150))
    print(f'end precision (processed locus): {end_precision(processed):.2f}')
    print(f'end precision (degradation):     {end_precision(decay):.2f}')
