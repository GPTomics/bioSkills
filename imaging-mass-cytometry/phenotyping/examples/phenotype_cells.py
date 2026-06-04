'''IMC phenotyping: cofactor-1 arcsinh, lineage-only clustering, double-positive diagnostic.

Encodes three IMC-specific decisions a CyTOF/scRNA pipeline gets wrong: the arcsinh
cofactor (1 not 5), clustering on lineage markers only (state markers fragment a type),
and treating an impossible double-positive as a segmentation/spillover artifact to be
diagnosed on the image, not a discovery.
'''
# Reference: scanpy 1.10+, anndata 0.10+, squidpy 1.3+, numpy 1.26+ | Verify API if version differs
import anndata as ad
import scanpy as sc
import squidpy as sq
import numpy as np

adata = ad.read_h5ad('imc_segmented.h5ad')
adata.layers['counts'] = adata.X.copy()
adata.X = np.arcsinh(adata.X / 1.0)   # cofactor ~1 for IMC means; 5 (CyTOF) over-compresses

# lineage markers only -- mixing state markers (Ki67, PD-1) splits one type into
# proliferating/resting pseudo-types
lineage = [m for m in ['CD45', 'CD3', 'CD8', 'CD4', 'CD20', 'CD68', 'E-cadherin'] if m in adata.var_names]
sub = adata[:, lineage].copy()
sc.pp.pca(sub, n_comps=min(15, len(lineage)))
sc.pp.neighbors(sub, n_neighbors=15)
sc.tl.leiden(sub, resolution=0.5)
adata.obs['leiden'] = sub.obs['leiden'].values
print(f'{adata.obs["leiden"].nunique()} clusters; validate with HELD-OUT evidence, not the clustering markers')

def _col(adata, marker):
    # dense 1-D values for one marker, robust to sparse X (np.asarray on a sparse matrix
    # returns a 0-d object array, not the values)
    x = adata[:, marker].X
    return (x.toarray() if hasattr(x, 'toarray') else np.asarray(x)).ravel()

def double_positive_is_artifact(adata, marker_a, marker_b, pct=75):
    # flag cells positive for a mutually-exclusive pair, then test whether they sit next to
    # a donor-type cell -- border/donor-adjacency means spillover or a merged segment
    a = _col(adata, marker_a)
    b = _col(adata, marker_b)
    suspect = (a > np.percentile(a, pct)) & (b > np.percentile(b, pct))
    sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True)
    conn = adata.obsp['spatial_connectivities']
    donor = b > np.percentile(b, pct)              # true B cells donate CD20
    neighbor_is_donor = (conn[suspect][:, donor].sum(axis=1) > 0)
    frac = float(np.asarray(neighbor_is_donor).mean()) if suspect.sum() else 0.0
    return suspect.sum(), frac

if 'CD3' in adata.var_names and 'CD20' in adata.var_names:
    n, frac = double_positive_is_artifact(adata, 'CD3', 'CD20')
    print(f'CD3+CD20+ suspects: {n}; fraction adjacent to a B cell: {frac:.2f} (high -> artifact)')

adata.write('imc_phenotyped.h5ad')
