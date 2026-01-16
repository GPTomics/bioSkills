'''Parse multiple trees from a single file'''

from Bio import Phylo
from io import StringIO

multi_tree_string = '''((A,B),(C,D));
((A,C),(B,D));
((A,D),(B,C));'''

trees = list(Phylo.parse(StringIO(multi_tree_string), 'newick'))
print(f'Loaded {len(trees)} trees\n')

for i, tree in enumerate(trees):
    terminals = [t.name for t in tree.get_terminals()]
    print(f'Tree {i+1}: {terminals}')
