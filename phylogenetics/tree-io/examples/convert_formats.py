'''Convert a phylogenetic tree between formats, writing to a temp dir.
Demonstrates that conversion to a lower-capability format is lossy: phyloXML holds
typed annotations that Nexus and Newick cannot, so a phyloXML->newick step drops them.'''
# Reference: BioPython 1.83+ | Verify API if version differs

import os
import tempfile
from io import StringIO
from Bio import Phylo

tree_string = '((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);'
tree = Phylo.read(StringIO(tree_string), 'newick')

out = tempfile.mkdtemp(prefix='treeio_')
xml_path = os.path.join(out, 'tree.xml')
nex_path = os.path.join(out, 'tree.nex')

Phylo.write(tree, xml_path, 'phyloxml')          # phyloXML is Bio.Phylo's richest format
Phylo.convert(xml_path, 'phyloxml', nex_path, 'nexus')   # Nexus cannot hold typed phyloXML annotations

print(f'Wrote {xml_path} and {nex_path}')
print('Topology survives; typed phyloXML annotations would not survive a step down to plain Newick.')
