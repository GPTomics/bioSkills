'''Empirical p-value via shuffling for pairwise alignment.

Reference: BioPython 1.83+ | Verify API if version differs
Reference: Pearson 2013 Curr Protoc Bioinf 3.1; Altschul-Erickson 1985 MBE.

Two shuffle strategies:
- Mononucleotide / mono-residue shuffle: preserves composition only
- Dinucleotide / di-residue shuffle: preserves composition + transition frequencies
  (Altschul-Erickson 1985 algorithm; better null for biased compositions)

For protein, mono-residue shuffle is the standard; dinucleotide is mainly
relevant for DNA where local context is meaningful (use the ushuffle library).
'''
import random
from Bio.Align import PairwiseAligner, substitution_matrices

def shuffle_seq(seq, preserve='mono'):
    '''Shuffle a sequence preserving the chosen statistic.'''
    if preserve == 'mono':
        chars = list(seq)
        random.shuffle(chars)
        return ''.join(chars)
    raise NotImplementedError('Use the ushuffle library for dinucleotide shuffling')

def empirical_pvalue(seq1, seq2, aligner, n_shuffles=1000, seed=42):
    '''Compute empirical p-value from a shuffled null distribution.

    Returns observed score, p-value, and the full null distribution.
    The (n_at_or_above + 1) / (n_shuffles + 1) correction (Phipson & Smyth 2010)
    avoids reporting p == 0 when no shuffle exceeds the observed score.
    '''
    random.seed(seed)
    observed = aligner.score(seq1, seq2)
    null_scores = [aligner.score(shuffle_seq(seq1), seq2) for _ in range(n_shuffles)]
    n_at_or_above = sum(1 for s in null_scores if s >= observed)
    pvalue = (n_at_or_above + 1) / (n_shuffles + 1)
    return observed, pvalue, null_scores

if __name__ == '__main__':
    aligner = PairwiseAligner(mode='local',
                              substitution_matrix=substitution_matrices.load('BLOSUM62'),
                              open_gap_score=-11, extend_gap_score=-1)
    observed, p, _ = empirical_pvalue('MKTIIALSYIFCLVFA', 'MKAIIVCSCLLVFFA', aligner, n_shuffles=1000)
    print(f'Observed: {observed:.1f}  Empirical p-value: {p:.4f}')
