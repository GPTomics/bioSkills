'''Global pairwise alignment of two DNA sequences'''

from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

seq1 = Seq('ACCGGTAACGTAG')
seq2 = Seq('ACCGTTAACGAAG')

aligner = PairwiseAligner(mode='global', match_score=2, mismatch_score=-1, open_gap_score=-10, extend_gap_score=-0.5)

alignments = aligner.align(seq1, seq2)
print(f'Found {len(alignments)} optimal alignment(s)')
print(f'Score: {alignments[0].score}\n')
print(alignments[0])
