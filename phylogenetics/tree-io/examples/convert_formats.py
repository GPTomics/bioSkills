'''Convert phylogenetic tree between formats'''

from Bio import Phylo
from io import StringIO

tree_string = '((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);'
tree = Phylo.read(StringIO(tree_string), 'newick')

Phylo.write(tree, 'output.xml', 'phyloxml')
print('Converted Newick to PhyloXML')

Phylo.convert('output.xml', 'phyloxml', 'output.nex', 'nexus')
print('Converted PhyloXML to Nexus')
