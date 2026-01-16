'''Build a Neighbor Joining tree from a multiple sequence alignment'''

from Bio import AlignIO, Phylo
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

calculator = DistanceCalculator('identity')
dm = calculator.get_distance(alignment)

print('Distance Matrix:')
print(dm)

constructor = DistanceTreeConstructor(calculator, 'nj')
tree = constructor.build_tree(alignment)
tree.ladderize()

print('\nNeighbor Joining Tree:')
Phylo.draw_ascii(tree)
