'''Enumerate, filter, and outcome-rank SpCas9 guide RNAs for gene knockout.

On-target activity scoring is intentionally NOT hand-rolled here: it is model- and
context-specific and belongs to a published model (route to CRISPOR, which selects
Rule Set 2/Azimuth for U6/lentiviral or CRISPRscan for T7/embryo per Haeussler 2016).
This script does the parts that ARE deterministic and local: PAM enumeration, hard
filters, and a microhomology-based out-of-frame estimate (Bae 2014 framework) that
predicts whether repair is frameshift-rich -- the property that actually drives
knockout. For a full repair-outcome distribution use inDelphi/FORECasT/Lindel.
'''
# Reference: biopython 1.83+ | Verify API if version differs

from Bio.Seq import Seq
import re
from math import exp

GUIDE_LEN = 20                 # SpCas9 spacer length
GC_MIN, GC_MAX = 0.40, 0.70    # outside this band on-target activity falls off (Doench 2014); soft penalty
MH_MIN = 3                     # shortest microhomology counted for MMEJ deletions
LENGTH_WEIGHT = 20.0           # deletion-length decay; longer MMEJ deletions are rarer (Bae 2014 framework)


def gc_fraction(spacer):
    return sum(c in 'GC' for c in spacer) / len(spacer)


def find_guides(sequence, guide_length=GUIDE_LEN):
    '''Enumerate SpCas9 (NGG) spacers on both strands; the spacer lies 5' of the PAM.

    cut is reported in forward-strand coordinates ~3 bp 5' of the PAM (Jinek 2012).
    '''
    seq = sequence.upper()
    n = len(seq)
    guides = []
    for m in re.finditer(r'(?=([ACGT]GG))', seq):
        pos = m.start()
        if pos >= guide_length:
            guides.append({'spacer': seq[pos - guide_length:pos], 'pam': seq[pos:pos + 3],
                           'cut': pos - 3, 'strand': '+'})
    rc = str(Seq(seq).reverse_complement())
    for m in re.finditer(r'(?=([ACGT]GG))', rc):
        pos = m.start()
        if pos >= guide_length:
            guides.append({'spacer': rc[pos - guide_length:pos], 'pam': rc[pos:pos + 3],
                           'cut': n - (pos - 3), 'strand': '-'})
    return guides


def passes_u6_filters(spacer):
    '''Hard filters for U6/H1 Pol-III expression. Skip these for in-vitro T7/RNP delivery.'''
    return 'TTTT' not in spacer and GC_MIN <= gc_fraction(spacer) <= GC_MAX   # TTTT terminates Pol III


def prepend_u6_g(spacer):
    '''U6 initiates best with a 5' G: prepend one (do NOT replace the first base).'''
    return spacer if spacer.startswith('G') else 'G' + spacer


def out_of_frame_score(local_seq, cut):
    '''Microhomology out-of-frame estimate (Bae 2014 framework).

    Enumerate pairs of identical microhomologies (>= MH_MIN) spanning the cut; each
    licenses an MMEJ deletion whose length is the distance between the copies. Weight
    each deletion by length decay and classify by frame. Returns the percent of weighted
    MMEJ patterns that are out-of-frame (frameshift); guides with >66 are preferred for
    reliable knockout. The EXACT published pattern weights live in the authors'
    implementation -- use it (or inDelphi/FORECasT) for production ranking; this is a
    transparent local approximation, not the canonical constants.
    '''
    left, right = local_seq[:cut], local_seq[cut:]
    weighted_oof, weighted_total = 0.0, 0.0
    seen = set()
    for k in range(min(len(left), len(right)), MH_MIN - 1, -1):
        for i in range(len(left) - k + 1):
            mh = left[i:i + k]
            for j in range(len(right) - k + 1):
                if right[j:j + k] != mh:
                    continue
                del_len = (cut + j) - i   # bases removed: gap between copies + one retained microhomology copy
                key = (i, cut + j)
                if del_len <= 0 or key in seen:
                    continue
                seen.add(key)
                gc = sum(c in 'GC' for c in mh)
                weight = exp(-del_len / LENGTH_WEIGHT) * ((k - gc) + 2 * gc)
                weighted_total += weight
                if del_len % 3 != 0:
                    weighted_oof += weight
    return 100.0 * weighted_oof / weighted_total if weighted_total else 0.0


def design_guides(sequence, delivery='u6', n_guides=5):
    '''Return candidate guides with filters applied and an out-of-frame estimate.

    delivery='u6' applies the Pol-III hard filters; delivery='t7'/'rnp' skips them.
    Final on-target ranking should come from CRISPOR (context-valid model), not from here.
    '''
    seq = sequence.upper()
    guides = find_guides(seq)
    out = []
    for g in guides:
        if delivery == 'u6' and not passes_u6_filters(g['spacer']):
            continue
        g['gc'] = gc_fraction(g['spacer'])
        g['expression_spacer'] = prepend_u6_g(g['spacer']) if delivery == 'u6' else g['spacer']
        lo, hi = max(0, g['cut'] - 30), min(len(seq), g['cut'] + 30)
        g['out_of_frame'] = out_of_frame_score(seq[lo:hi], g['cut'] - lo)
        out.append(g)
    out.sort(key=lambda x: x['out_of_frame'], reverse=True)
    return out[:n_guides]


if __name__ == '__main__':
    target = ('ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGT'
              'GTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTG'
              'CATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGATATAACCAAA')

    guides = design_guides(target, delivery='u6', n_guides=5)
    print(f'{len(guides)} candidate guides (U6 filters applied), ranked by out-of-frame estimate:')
    for i, g in enumerate(guides, 1):
        print(f"{i}. {g['spacer']} PAM={g['pam']} strand={g['strand']} cut={g['cut']} "
              f"GC={g['gc']:.0%} out-of-frame~{g['out_of_frame']:.0f} "
              f"(rank on-target with CRISPOR; design 3-6 and validate)")
