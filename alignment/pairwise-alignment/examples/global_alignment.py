'''Global pairwise alignment examples with different scoring schemes'''

from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq

# DNA alignment with custom scoring
dna_seq1 = Seq('ACCGGTAACGTAG')
dna_seq2 = Seq('ACCGTTAACGAAG')

aligner = PairwiseAligner(mode='global', match_score=2, mismatch_score=-1, open_gap_score=-10, extend_gap_score=-0.5)
alignments = aligner.align(dna_seq1, dna_seq2)
print(f'DNA alignment score: {alignments[0].score}')
print(alignments[0])

# Protein alignment - BLOSUM62 (default, moderate divergence)
prot_seq1 = Seq('MKFLILLFNILCLFPVLAADNH')
prot_seq2 = Seq('MKFLVLLFNILCLFPVLAADHH')

aligner = PairwiseAligner(mode='global')
aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
aligner.open_gap_score = -11
aligner.extend_gap_score = -1
alignments = aligner.align(prot_seq1, prot_seq2)
print(f'BLOSUM62 alignment score: {alignments[0].score}')
print(alignments[0])

# BLOSUM45 for distant/divergent sequences (more permissive)
aligner.substitution_matrix = substitution_matrices.load('BLOSUM45')
aligner.open_gap_score = -12
aligner.extend_gap_score = -2
alignments = aligner.align(prot_seq1, prot_seq2)
print(f'BLOSUM45 alignment score: {alignments[0].score}')

# BLOSUM80 for closely related sequences (more stringent)
aligner.substitution_matrix = substitution_matrices.load('BLOSUM80')
aligner.open_gap_score = -10
aligner.extend_gap_score = -1
alignments = aligner.align(prot_seq1, prot_seq2)
print(f'BLOSUM80 alignment score: {alignments[0].score}')

# Affine vs linear gap penalties comparison
aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
aligner.open_gap_score = -11
aligner.extend_gap_score = -1
affine_alignments = aligner.align(prot_seq1, prot_seq2)

aligner.open_gap_score = -5
aligner.extend_gap_score = -5
linear_alignments = aligner.align(prot_seq1, prot_seq2)
print(f'Affine gap score: {affine_alignments[0].score}, Linear gap score: {linear_alignments[0].score}')
