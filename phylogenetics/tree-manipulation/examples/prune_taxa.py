'''Prune taxa while preserving every surviving patristic distance.

Dropping a tip leaves a degree-2 'knee' node; the correct behavior suppresses it
and ADDS its branch length to the child so path lengths are conserved. Bio.Phylo
prune does this by default. The script proves it by checking a pairwise distance
before and after pruning.
'''
# Reference: biopython 1.83+ | Verify API if version differs

from Bio import Phylo
from io import StringIO

tree_string = '((Human:0.1,Chimp:0.2):0.3,(Mouse:0.4,Rat:0.5):0.6,Zebrafish:1.0);'

tree = Phylo.read(StringIO(tree_string), 'newick')
print('Original tree:')
Phylo.draw_ascii(tree)

dist_before = tree.distance('Human', 'Mouse')

keep = {'Human', 'Mouse', 'Rat'}
for term in list(tree.get_terminals()):
    if term.name not in keep:
        tree.prune(term)

dist_after = tree.distance('Human', 'Mouse')
print('\nAfter keeping only Human, Mouse, Rat:')
Phylo.draw_ascii(tree)
print(f'Human-Mouse distance: before={dist_before}, after={dist_after} (must be unchanged)')
