'''Calculate pairwise distances between taxa in a tree'''

from Bio import Phylo
from io import StringIO

tree_string = '((Human:0.1,Chimp:0.2):0.3,(Mouse:0.4,Rat:0.5):0.6);'
tree = Phylo.read(StringIO(tree_string), 'newick')

print('Tree:')
Phylo.draw_ascii(tree)

print('\nPairwise distances (sum of branch lengths):')
terminals = tree.get_terminals()
for i, t1 in enumerate(terminals):
    for t2 in terminals[i+1:]:
        d = tree.distance(t1, t2)
        print(f'{t1.name} - {t2.name}: {d:.2f}')

print(f'\nTotal branch length: {tree.total_branch_length():.2f}')

depths = tree.depths()
print(f'Tree height: {max(depths.values()):.2f}')
