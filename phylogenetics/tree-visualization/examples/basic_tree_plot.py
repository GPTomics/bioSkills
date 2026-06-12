'''Draw a phylogram and export to a vector PDF in a temp dir (no CWD strays)'''
# Reference: biopython 1.83+, matplotlib 3.8+ | Verify API if version differs

import os
import tempfile
from io import StringIO
from Bio import Phylo
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

tree_string = '((Human:0.1,Chimp:0.2):0.3,(Mouse:0.4,Rat:0.5):0.6);'
tree = Phylo.read(StringIO(tree_string), 'newick')
tree.ladderize()                                  # legibility only; ordering carries NO phylogenetic meaning

fig, ax = plt.subplots(figsize=(10, 6))
Phylo.draw(tree, axes=ax, do_show=False)
ax.set_title('Example phylogram (branch length = subs/site)')

out_dir = tempfile.mkdtemp(prefix='tree_viz_')
out_path = os.path.join(out_dir, 'basic_tree.pdf')   # vector keeps text/lines sharp at print size
fig.savefig(out_path, bbox_inches='tight')
plt.close(fig)
print('Saved to', out_path)
