'''Read and inspect a Newick format phylogenetic tree'''

from Bio import Phylo
from io import StringIO

tree_string = '((Human:0.1,Chimp:0.2):0.3,(Mouse:0.4,Rat:0.5):0.6);'
tree = Phylo.read(StringIO(tree_string), 'newick')

print('Tree structure:')
Phylo.draw_ascii(tree)

print(f'\nTerminal taxa: {[t.name for t in tree.get_terminals()]}')
print(f'Total branch length: {tree.total_branch_length():.2f}')
print(f'Is bifurcating: {tree.is_bifurcating()}')
