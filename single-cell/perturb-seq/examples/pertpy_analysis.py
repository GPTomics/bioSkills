'''Perturb-seq analysis with Pertpy'''
import scanpy as sc
import pertpy as pt
import pandas as pd

# Load data
adata = sc.read_h5ad('perturb_seq.h5ad')

# Standard preprocessing
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)

# Assign perturbations from guide calls
# Assumes guide_id column exists in obs
# Non-targeting guides labeled as 'NT' or 'non-targeting'
print(f'Perturbations: {adata.obs["guide_id"].nunique()}')
print(f'Non-targeting cells: {(adata.obs["guide_id"] == "NT").sum()}')

# Differential expression per perturbation vs non-targeting
de = pt.tl.PseudobulkDE(adata)
de.fit(
    groupby='guide_id',
    control='NT',  # Non-targeting control
    layer=None  # Use .X
)

results = de.results()
# FDR < 0.05, |log2FC| > 1
sig_results = results[(results['pvals_adj'] < 0.05) & (abs(results['logfoldchanges']) > 1)]
print(f'Significant DE genes across perturbations: {len(sig_results)}')

# Compute perturbation signatures
# Embedding distance from control cells
ps = pt.tl.PerturbationSignature(adata)
ps.compute(
    groupby='guide_id',
    control='NT'
)

# Get signature scores per perturbation
signatures = ps.get_signatures()
signatures.to_csv('perturbation_signatures.csv')

# Cluster perturbations by effect
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, color='guide_id', save='_perturbations.png')
