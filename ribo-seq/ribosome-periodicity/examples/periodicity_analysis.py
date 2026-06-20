'''Score 3-nucleotide periodicity from CDS-body P-site coverage.

The defensible metric is the frame-0 fraction over the CDS body. The start
(initiation) and stop (termination) peaks dwarf the body and carry their own
phase, so they are trimmed before scoring. A spectral check (power at period 3)
is reported as a secondary signal. For real P-site offset calibration use
riboWaltz (psite/psite_info) or plastid's metagene+psite CLI scripts; this
script demonstrates the scoring once a per-nucleotide P-site coverage vector exists.
'''
# Reference: numpy 1.26+, scipy 1.12+ | Verify API if version differs

import numpy as np
from scipy.fft import rfft, rfftfreq

# Trim windows: start-codon initiation peak spans the first ~15 codons; the stop
# (termination) peak the last ~5 codons. Trimming in nucleotides (codons x 3).
TRIM_START_NT = 45
TRIM_STOP_NT = 15
# Frame-0 fraction thresholds (fraction of in-CDS body P-sites in frame 0; null = 0.33)
GOOD_FRAME0 = 0.60
MARGINAL_FRAME0 = 0.45


def body_frame_fraction(psite_coverage, trim_start=TRIM_START_NT, trim_stop=TRIM_STOP_NT):
    '''Fraction of CDS-body P-sites in each reading frame; frame 0 should dominate.'''
    body = psite_coverage[trim_start:len(psite_coverage) - trim_stop]
    frames = np.array([body[f::3].sum() for f in range(3)], dtype=float)
    total = frames.sum()
    return frames / total if total else frames


def period3_power(psite_coverage, trim_start=TRIM_START_NT, trim_stop=TRIM_STOP_NT):
    '''Relative spectral power at period 3 over CDS-body coverage (secondary metric).

    Runs on uniform body coverage, NOT the start-codon metagene (which is one peak).
    '''
    body = psite_coverage[trim_start:len(psite_coverage) - trim_stop]
    body = body - body.mean()
    if not np.any(body):
        return 0.0
    power = np.abs(rfft(body)) ** 2
    freq = rfftfreq(len(body))
    idx3 = np.argmin(np.abs(freq - 1.0 / 3.0))
    return float(power[idx3] / power.sum())


def classify(frame0_fraction):
    if frame0_fraction >= GOOD_FRAME0:
        return 'good'
    if frame0_fraction >= MARGINAL_FRAME0:
        return 'marginal'
    return 'poor'


def _simulate_cds(n_codons=300, frame0_weight=0.75, depth=40, seed=0):
    '''Simulate per-nt P-site coverage with frame-0 enrichment plus start/stop peaks.'''
    rng = np.random.default_rng(seed)
    n_nt = n_codons * 3
    base = np.zeros(n_nt)
    in_frame = rng.poisson(depth * frame0_weight, n_codons)
    off_frame = rng.poisson(depth * (1 - frame0_weight) / 2, (n_codons, 2))
    base[0::3][:n_codons] = in_frame
    base[1::3][:n_codons] = off_frame[:, 0]
    base[2::3][:n_codons] = off_frame[:, 1]
    base[0:3] += depth * 30      # initiation peak
    base[-6:-3] += depth * 12    # termination peak
    return base


if __name__ == '__main__':
    print('Ribosome periodicity scoring (simulated CDS-body coverage)')
    for label, weight in [('good library', 0.78), ('marginal', 0.50), ('no periodicity', 0.34)]:
        cov = _simulate_cds(frame0_weight=weight, seed=hash(label) % 100)
        frames = body_frame_fraction(cov)
        f0 = frames[0]
        print(f'  {label:16s} frame fractions={np.round(frames, 3)} '
              f'frame0={f0:.3f} ({classify(f0)}) period3_power={period3_power(cov):.3f}')
