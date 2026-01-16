'''Draw a basic phylogenetic tree and save to file'''

from Bio import Phylo
from io import StringIO
import matplotlib.pyplot as plt

tree_string = '((Human:0.1,Chimp:0.2):0.3,(Mouse:0.4,Rat:0.5):0.6);'
tree = Phylo.read(StringIO(tree_string), 'newick')
tree.ladderize()

fig, ax = plt.subplots(figsize=(10, 6))
Phylo.draw(tree, axes=ax, do_show=False)
ax.set_title('Example Phylogenetic Tree')
plt.savefig('basic_tree.png', dpi=300, bbox_inches='tight')
plt.close()
print('Saved to basic_tree.png')
