'''Position-specific scoring matrix (PSSM) with Laplace pseudocounts.

Reference: BioPython 1.83+ | Verify API if version differs
Reference: Henikoff JG & Henikoff S 1996 Bioinf 12:135-143 (data-dependent pseudocounts).

Without pseudocounts, log-odds against background diverge to negative infinity
at any column missing a residue. Laplace add-one is the minimal correct
smoothing for production use; HMMER uses Dirichlet mixtures for finer control.

Default background is Robinson & Robinson 1991 (protein only); A/C/G/T are present
in both protein and DNA alphabets, so this script will silently run on DNA input.
For DNA alignments, pass background={'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
explicitly.
'''
import math
from collections import Counter
from Bio import AlignIO

ROBINSON_BACKGROUND = {
    'A': 0.0780, 'R': 0.0512, 'N': 0.0427, 'D': 0.0530, 'C': 0.0193,
    'Q': 0.0419, 'E': 0.0629, 'G': 0.0738, 'H': 0.0224, 'I': 0.0526,
    'L': 0.0922, 'K': 0.0596, 'M': 0.0224, 'F': 0.0399, 'P': 0.0508,
    'S': 0.0712, 'T': 0.0584, 'W': 0.0133, 'Y': 0.0327, 'V': 0.0653,
}

def pssm_with_pseudocounts(alignment, background=None, pseudocount=1.0):
    if background is None:
        background = ROBINSON_BACKGROUND
    pssm = []
    for col_idx in range(alignment.get_alignment_length()):
        column = alignment[:, col_idx].replace('-', '')
        n = len(column)
        counts = Counter(column)
        pssm.append({
            residue: math.log2(((counts.get(residue, 0) + pseudocount * background[residue]) / (n + pseudocount))
                                / background[residue])
            for residue in background
        })
    return pssm

def score_site(sequence, pssm):
    return sum(pssm[i].get(residue, 0) for i, residue in enumerate(sequence))

if __name__ == '__main__':
    alignment = AlignIO.read('alignment.fasta', 'fasta')
    pssm = pssm_with_pseudocounts(alignment)
    print(f'PSSM with {len(pssm)} positions')
    for i, col_scores in enumerate(pssm[:5]):
        top = sorted(col_scores.items(), key=lambda kv: -kv[1])[:3]
        print(f'  Position {i}: {", ".join(f"{r}={s:+.2f}" for r, s in top)}')
