'''Collapse low-support branches into SOFT (uncertainty) polytomies.

Bio.Phylo parses an internal-node support label into clade.confidence. Collapsing
branches below a support cutoff produces soft polytomies meaning 'order unresolved'
-- NOT hard polytomies meaning 'simultaneous radiation'. The cutoff depends on the
support scale: 70 for standard bootstrap, 95 for UFBoot2 (different scales).
'''
# Reference: biopython 1.83+ | Verify API if version differs

from Bio import Phylo
from io import StringIO

# Internal nodes carry bootstrap support after the closing paren: (A,B)40 etc.
tree_string = '(((Human:0.1,Chimp:0.1)40:0.2,Gorilla:0.3)95:0.4,Orangutan:0.5);'

tree = Phylo.read(StringIO(tree_string), 'newick')
print('Internal-node support values:')
for clade in tree.get_nonterminals():
    print(f'  confidence={clade.confidence}')

bootstrap_cutoff = 70   # use 95 for UFBoot2; these are not the same scale
tree.collapse_all(lambda c: c.confidence is not None and c.confidence < bootstrap_cutoff)

print(f'\nAfter collapsing branches with bootstrap < {bootstrap_cutoff} (soft polytomy = unresolved):')
Phylo.draw_ascii(tree)
