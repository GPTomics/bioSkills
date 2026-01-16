'''Protein alignment using BLOSUM62 substitution matrix'''

from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq

protein1 = Seq('MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQQIAAALEHHHHHH')
protein2 = Seq('MKTAYIAKQRQISFVKSHFSRQLEERLDLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQQIA')

blosum62 = substitution_matrices.load('BLOSUM62')
aligner = PairwiseAligner(mode='global', substitution_matrix=blosum62, open_gap_score=-11, extend_gap_score=-1)

alignments = aligner.align(protein1, protein2)
print(f'Score: {alignments[0].score}\n')
print(alignments[0])
