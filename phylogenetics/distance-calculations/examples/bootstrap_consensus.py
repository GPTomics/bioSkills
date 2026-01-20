'''Generate bootstrap consensus tree from alignment'''

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo.Consensus import bootstrap_consensus, majority_consensus
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sequences = [
    SeqRecord(Seq('ATGCATGCATGCATGC'), id='Human'),
    SeqRecord(Seq('ATGCATGCATGAATGC'), id='Chimp'),
    SeqRecord(Seq('ATGCATGAATGCATGC'), id='Gorilla'),
    SeqRecord(Seq('ATGAATGCATGCATGC'), id='Mouse'),
    SeqRecord(Seq('ATGAATGAATGCATGC'), id='Rat'),
]
alignment = MultipleSeqAlignment(sequences)

calculator = DistanceCalculator('identity')
constructor = DistanceTreeConstructor(calculator, 'nj')

print('Building bootstrap consensus (50 replicates)...')
consensus_tree = bootstrap_consensus(alignment, 50, constructor, majority_consensus)
consensus_tree.ladderize()

print('\nMajority Rule Consensus Tree:')
Phylo.draw_ascii(consensus_tree)

# Expected output: Bootstrap values appear as branch confidences (0-100%)
# Interpretation: <50 = weak, 50-70 = moderate, 70-90 = good, >90 = strong support
# For publication: typically require >70% bootstrap support to report a clade
