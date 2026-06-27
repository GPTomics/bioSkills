'''Lineage tree reconstruction with Cassiopeia'''
# Reference: cassiopeia 2.0+, matplotlib 3.8+, numpy 1.26+, scanpy 1.10+ | Verify API if version differs
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

print(f'cells {tree.n_cell}  characters {tree.n_character}  missing {(char_matrix == -1).mean():.2%}')

# HybridSolver is the scalable default: greedy top split, ILP on small subclades
# collapse_mutationless_edges removes internal edges with no supporting mutation
solver = cas.solver.HybridSolver(
    top_solver=cas.solver.VanillaGreedySolver(),
    bottom_solver=cas.solver.ILPSolver(),
    cell_cutoff=200
)
solver.solve(tree, collapse_mutationless_edges=True)

newick = tree.get_newick()

# Infer ancestral states on internal nodes
tree.reconstruct_ancestral_characters()

# Compare against a neighbor-joining tree to gauge topology robustness
nj_tree = cas.data.CassiopeiaTree(character_matrix=char_matrix, cell_meta=cell_meta)
cas.solver.NeighborJoiningSolver(
    dissimilarity_function=cas.solver.dissimilarity_functions.weighted_hamming_distance
).solve(nj_tree)
rf, rf_max = cas.critique.robinson_foulds(tree, nj_tree)
print(f'Robinson-Foulds {rf}/{rf_max}  triplets-correct {cas.critique.triplets_correct(tree, nj_tree)}')

if 'cell_type' in cell_meta.columns:
    cas.pl.plot_matplotlib(tree, meta_data=['cell_type'])
