import snapatac2 as snap
import scanpy as sc
import numpy as np

# Read 10X scATAC-seq data
# fragment file is the primary input for scATAC analysis
adata = snap.pp.import_data(
    'fragments.tsv.gz',
    chrom_sizes=snap.genome.hg38,
    file='scatac_processed.h5ad',
    min_num_fragments=1000,
    sorted_by_barcode=False
)

# Add cell metadata if available
# metadata = pd.read_csv('singlecell.csv', index_col=0)
# adata.obs = adata.obs.join(metadata)

# QC metrics
snap.metrics.tsse(adata, snap.genome.hg38)  # TSS enrichment score
snap.metrics.frag_size_distr(adata)  # Fragment size distribution

# QC filtering
# min_tsse > 2: minimum TSS enrichment (chromatin accessibility signal)
# min_counts > 3000: minimum fragments per cell
# max_counts < 50000: remove potential doublets or debris
adata = adata[(adata.obs['tsse'] > 2) & (adata.obs['n_fragment'] > 3000) & (adata.obs['n_fragment'] < 50000)].copy()
print(f'Cells after QC: {adata.n_obs}')

# Create tile matrix (500bp bins)
snap.pp.add_tile_matrix(adata, bin_size=500)

# Feature selection - select top variable features
snap.pp.select_features(adata, n_features=50000)

# Dimensionality reduction with spectral embedding (similar to LSI in Signac)
snap.tl.spectral(adata, n_comps=30)

# Check depth correlation with first component
# If correlated, exclude from downstream analysis (similar to Signac)
depth_cor = np.corrcoef(adata.obs['n_fragment'], adata.obsm['X_spectral'][:, 0])[0, 1]
use_dims = range(1, 30) if abs(depth_cor) > 0.5 else range(30)

# UMAP and clustering
snap.tl.umap(adata, use_dims=use_dims)
snap.pp.knn(adata, use_dims=use_dims)
snap.tl.leiden(adata, resolution=0.5)

# Save UMAP plot
sc.settings.figdir = './'
sc.pl.umap(adata, color='leiden', save='_scatac.pdf')

# Gene activity matrix (gene accessibility scores)
gene_matrix = snap.pp.make_gene_matrix(adata, snap.genome.hg38)
gene_matrix.var_names_make_unique()

# Normalize gene activity for visualization
sc.pp.normalize_total(gene_matrix, target_sum=1e4)
sc.pp.log1p(gene_matrix)

# Plot marker gene activities
marker_genes = ['CD34', 'MS4A1', 'CD3D', 'CD14', 'GATA1', 'PAX5']
available_markers = [g for g in marker_genes if g in gene_matrix.var_names]
if available_markers:
    sc.pl.umap(gene_matrix, color=available_markers, save='_gene_activity.pdf', ncols=3)

# Save processed data
adata.write('scatac_processed.h5ad')
print('Done. Saved to scatac_processed.h5ad')
