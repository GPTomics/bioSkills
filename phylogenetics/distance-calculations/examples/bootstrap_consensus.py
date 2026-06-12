'''Generate bootstrap consensus tree from alignment'''
# Reference: biopython 1.83+ | Verify API if version differs

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

calculator = DistanceCalculator('identity')   # identity-only p-distance; ape dist.dna needed for model correction
constructor = DistanceTreeConstructor(calculator, 'nj')

print('Building bootstrap consensus (50 replicates)...')
consensus_tree = bootstrap_consensus(alignment, 50, constructor, majority_consensus)
consensus_tree.ladderize()

print('\nMajority Rule Consensus Tree:')
Phylo.draw_ascii(consensus_tree)

# Bootstrap values appear as branch confidences (0-100%): <50 weak, 50-70 moderate, 70-90 good, >90 strong.
# Caveat: support is sampling PRECISION, not accuracy -- a bias in the distances reproduces every replicate.
