'''Draw a tree with tip labels and support labeled by its measure, export to temp SVG'''
# Reference: biopython 1.83+, matplotlib 3.8+ | Verify API if version differs

import os
import tempfile
from io import StringIO
from Bio import Phylo
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

tree_string = '((A:0.15,B:0.22)90:0.08,(C:0.35,D:0.41)72:0.12);'   # internal labels are bootstrap percentages
tree = Phylo.read(StringIO(tree_string), 'newick')
tree.ladderize()

def tip_only(clade):
    return clade.name if clade.is_terminal() else ''

def support_label(clade):
    # the measure MUST be named in the caption; bare integers default-read as bootstrap and over-read other scales
    if not clade.is_terminal() and clade.confidence is not None:
        return f'{clade.confidence:.0f}'
    return ''

fig, ax = plt.subplots(figsize=(10, 6))
Phylo.draw(tree, axes=ax, do_show=False, label_func=tip_only, branch_labels=support_label)
ax.set_title('Bootstrap support shown at internal nodes')

out_dir = tempfile.mkdtemp(prefix='tree_viz_')
out_path = os.path.join(out_dir, 'labeled_tree.svg')   # vector for publication
fig.savefig(out_path, bbox_inches='tight')
plt.close(fig)
print('Saved to', out_path)
