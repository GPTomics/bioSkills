'''Lineage tree reconstruction with Cassiopeia'''
import cassiopeia as cas
import pandas as pd
import numpy as np

# Option 1: Load pre-built character matrix
# Rows = cells, Columns = barcode sites
# Values: 0 = unedited, 1-N = different mutations, -1 = missing
char_matrix = pd.read_csv('character_matrix.csv', index_col=0)
cell_meta = pd.read_csv('cell_metadata.csv', index_col=0)

# Create CassiopeiaTree object
tree = cas.data.CassiopeiaTree(
    character_matrix=char_matrix,
    cell_meta=cell_meta
)

print(f'Cells: {tree.n_cell}')
print(f'Characters: {tree.n_character}')

# Tree reconstruction with Greedy solver
# Fast, suitable for most datasets
solver = cas.solver.VanillaGreedySolver()
solver.solve(tree)

# Alternative: Neighbor-joining (faster, less accurate)
# solver = cas.solver.NeighborJoiningSolver()

# Alternative: ILP solver (optimal but slow for large trees)
# solver = cas.solver.ILPSolver()

# Get tree in Newick format
newick = tree.get_newick()
with open('lineage_tree.newick', 'w') as f:
    f.write(newick)

# Compute tree statistics
# Tree depth, balance, etc.
tree.compute_tree_statistics()
print(f'Tree depth: {tree.tree_statistics["depth"]}')

# Annotate internal nodes
# Infer ancestral states
tree.reconstruct_ancestral_characters()

# Visualize with cell type annotations
if 'cell_type' in cell_meta.columns:
    cas.pl.tree_plot(
        tree,
        color_column='cell_type',
        save='lineage_tree_celltypes.png'
    )
