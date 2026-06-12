'''Extract a clade vs an induced subtree, and the non-monophyly trap.

common_ancestor returns the MRCA even when the targets are NOT monophyletic, in
which case the MRCA clade contains EXTRA taxa beyond the ones requested. Extracting
'the clade of X,Y,Z' is only valid when X,Y,Z are monophyletic; otherwise prune to
the taxon set to get the induced subtree instead.
'''
# Reference: biopython 1.83+ | Verify API if version differs

from Bio import Phylo
from io import StringIO

tree_string = '(((Human:0.1,Chimp:0.2):0.1,Gorilla:0.3):0.2,(Mouse:0.4,Rat:0.5):0.6);'

tree = Phylo.read(StringIO(tree_string), 'newick')
print('Full tree:')
Phylo.draw_ascii(tree)

primates = [tree.find_any(name='Human'), tree.find_any(name='Chimp'), tree.find_any(name='Gorilla')]
print(f'\nAre Human, Chimp, Gorilla monophyletic? {tree.is_monophyletic(primates)}')

clade = tree.common_ancestor({'name': 'Human'}, {'name': 'Gorilla'})
print(f'MRCA clade taxa: {[t.name for t in clade.get_terminals()]}')

# Non-monophyletic request: MRCA clade pulls in extra taxa
mixed = tree.common_ancestor({'name': 'Human'}, {'name': 'Mouse'})
print(f'\nMRCA of Human + Mouse pulls in extra taxa: {[t.name for t in mixed.get_terminals()]}')
print('For a non-monophyletic set, prune to the taxon set for the induced subtree instead.')
