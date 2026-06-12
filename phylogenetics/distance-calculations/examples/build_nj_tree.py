'''Build a Neighbor Joining tree from an alignment. Deliverable: a fast NJ tree -- but the
distance here is IDENTITY-ONLY (a p-distance), NOT a model correction. For divergent DNA use
ape dist.dna(model='TN93') for a model-corrected matrix; the matrix is the limitation, not NJ.'''
# Reference: biopython 1.83+ | Verify API if version differs

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sequences = [
    SeqRecord(Seq('ATGCATGCATGC'), id='Human'),
    SeqRecord(Seq('ATGCATGCATGA'), id='Chimp'),
    SeqRecord(Seq('ATGCATGAATGC'), id='Gorilla'),
    SeqRecord(Seq('ATGAATGCATGC'), id='Mouse'),
    SeqRecord(Seq('ATGAATGAATGC'), id='Rat'),
]
alignment = MultipleSeqAlignment(sequences)

calculator = DistanceCalculator('identity')   # identity-only: p-distance, NO multiple-hit correction
dm = calculator.get_distance(alignment)        # divergent data needs ape dist.dna for JC/K80/TN93/LogDet

print('Distance Matrix (uncorrected p-distance):')
print(dm)

constructor = DistanceTreeConstructor(calculator, 'nj')
tree = constructor.build_tree(alignment)       # NJ is consistent ONLY if the distances are correct/additive
tree.ladderize()

print('\nNeighbor Joining Tree:')
Phylo.draw_ascii(tree)
