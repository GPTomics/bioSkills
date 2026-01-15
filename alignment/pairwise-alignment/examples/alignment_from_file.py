'''Align sequences from a FASTA file'''

from Bio import SeqIO
from Bio.Align import PairwiseAligner

fasta_file = 'sequences.fasta'
records = list(SeqIO.parse(fasta_file, 'fasta'))

if len(records) < 2:
    print('Need at least 2 sequences to align')
else:
    seq1, seq2 = records[0].seq, records[1].seq

    aligner = PairwiseAligner(mode='global', match_score=2, mismatch_score=-1, open_gap_score=-10, extend_gap_score=-0.5)
    alignments = aligner.align(seq1, seq2)

    print(f'Aligning {records[0].id} vs {records[1].id}')
    print(f'Score: {alignments[0].score}\n')
    print(alignments[0])
