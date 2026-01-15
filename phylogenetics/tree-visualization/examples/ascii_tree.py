'''Display ASCII representation of a tree'''

from Bio import Phylo
from io import StringIO

tree_string = '((Human:0.1,Chimp:0.2):0.3,(Mouse:0.4,Rat:0.5):0.6,Zebrafish:1.0);'
tree = Phylo.read(StringIO(tree_string), 'newick')

print('Tree structure:')
print(tree)
print()
print('ASCII diagram:')
Phylo.draw_ascii(tree)
