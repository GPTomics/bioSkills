'''Local alignment to find best matching regions'''

from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

long_seq = Seq('NNNNNACCGGTAACGTAGNNNNNNNN')
short_seq = Seq('ACCGGTAACGTAG')

aligner = PairwiseAligner(mode='local', match_score=2, mismatch_score=-1, open_gap_score=-10, extend_gap_score=-0.5)

alignments = aligner.align(long_seq, short_seq)
print(f'Found {len(alignments)} optimal alignment(s)')
print(f'Score: {alignments[0].score}\n')
print(alignments[0])

# Show aligned coordinates
alignment = alignments[0]
print(f'\nAligned regions:')
print(f'Target: positions {alignment.aligned[0]}')
print(f'Query: positions {alignment.aligned[1]}')
