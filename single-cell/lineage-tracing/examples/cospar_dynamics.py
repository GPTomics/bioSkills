'''Clonal dynamics analysis with CoSpar'''
import cospar as cs
import scanpy as sc

# Load AnnData with lineage information
# Requires clone_id or barcode in obs
adata = sc.read_h5ad('lineage_traced.h5ad')

# Standard preprocessing
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# Initialize CoSpar
cs.settings.set_figure_params()

# Prepare clonal data
# Converts clone IDs to sparse matrix format
cs.tl.prepare_clones(
    adata,
    clone_key='clone_id',  # Column with barcode/clone ID
    min_clone_size=2       # Min cells per clone
)

# Infer transition map (fate probabilities)
# smooth_array: smoothing parameters for refinement
cs.tl.infer_Tmap(
    adata,
    smooth_array=[15, 10, 5],  # Multi-scale smoothing
    max_iter=3
)

# Compute fate map from progenitor to terminal state
# Requires cell_type annotation
cs.tl.fate_map(
    adata,
    source='Progenitor',  # Starting cell type
    map_key='Tmap'
)

# Visualize fate probabilities
cs.pl.fate_map(adata, source='Progenitor')

# Clone expansion analysis
# Track how individual clones grow over time
if 'time_point' in adata.obs.columns:
    cs.tl.clone_dynamics(
        adata,
        time_key='time_point',
        clone_key='clone_id'
    )
    cs.pl.clone_dynamics(adata)

# Export results
adata.write('cospar_results.h5ad')
