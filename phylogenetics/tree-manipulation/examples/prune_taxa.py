'''Prune taxa from a phylogenetic tree'''

from Bio import Phylo
from io import StringIO

tree_string = '((Human:0.1,Chimp:0.2):0.3,(Mouse:0.4,Rat:0.5):0.6,Zebrafish:1.0);'
tree = Phylo.read(StringIO(tree_string), 'newick')

print('Original tree:')
Phylo.draw_ascii(tree)
print(f'Taxa: {[t.name for t in tree.get_terminals()]}')

keep_taxa = {'Human', 'Chimp', 'Mouse'}
for term in tree.get_terminals():
    if term.name not in keep_taxa:
        tree.prune(term)

print('\nAfter keeping only Human, Chimp, Mouse:')
Phylo.draw_ascii(tree)
print(f'Taxa: {[t.name for t in tree.get_terminals()]}')
