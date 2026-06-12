'''Root a tree, treating the root as a separate inference, not a display choice.

Outgroup rooting (multiple close taxa) is preferred and checks ingroup monophyly.
Midpoint rooting assumes a clock and is a clock-limited fallback. For deep trees
with rate variation, prefer outgroup-free MAD / MinVar CLIs (run on the Newick file).
'''
# Reference: biopython 1.83+ | Verify API if version differs

from Bio import Phylo
from io import StringIO

tree_string = '(((Human:0.1,Chimp:0.2):0.3,Gorilla:0.4):0.2,(OutA:0.5,OutB:0.6):0.7);'

tree = Phylo.read(StringIO(tree_string), 'newick')
print('Unrooted/arbitrarily-rooted input:')
Phylo.draw_ascii(tree)

outgroup_clades = [tree.find_any(name='OutA'), tree.find_any(name='OutB')]
if tree.is_monophyletic(outgroup_clades):
    tree.root_with_outgroup({'name': 'OutA'}, {'name': 'OutB'})
    print('\nRooted with the monophyletic outgroup (OutA, OutB):')
    Phylo.draw_ascii(tree)
else:
    print('\nOutgroup is not monophyletic: root placement is unreliable')

midpoint_tree = Phylo.read(StringIO(tree_string), 'newick')
midpoint_tree.root_at_midpoint()
print('\nMidpoint-rooted (clock assumption; a long branch can slide the root):')
Phylo.draw_ascii(midpoint_tree)
