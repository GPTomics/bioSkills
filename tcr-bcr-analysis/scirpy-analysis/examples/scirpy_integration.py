'''Single-cell TCR/BCR analysis with scirpy'''

import scirpy as ir
import scanpy as sc
import pandas as pd

# Load scRNA-seq data
adata = sc.read_h5ad('scrnaseq.h5ad')
print(f'Loaded {adata.n_obs} cells, {adata.n_vars} genes')

# Add 10x VDJ data
# filtered_contig_annotations.csv from Cell Ranger VDJ
ir.io.read_10x_vdj(adata, 'filtered_contig_annotations.csv')
print(f'Cells with TCR/BCR: {adata.obs["has_ir"].sum()}')

# Chain quality control
# Identifies potential issues: doublets, orphan chains, ambiguous pairing
ir.tl.chain_qc(adata)

# QC categories:
# - multichain: >2 chains detected (doublet marker)
# - orphan: Only alpha OR beta, not both
# - extra VDJ/VJ: Additional chains beyond the pair
print('\nChain pairing QC:')
print(adata.obs['chain_pairing'].value_counts())

# Filter to high-quality cells
# Keep cells with proper chain pairing (single pair)
adata_clean = adata[adata.obs['chain_pairing'].isin(['single pair'])].copy()
print(f'\nRetained {adata_clean.n_obs} cells with proper pairing')

# Define clonotypes
# Uses amino acid CDR3 sequence identity
# cutoff=0 means exact match required
ir.pp.ir_dist(adata_clean, metric='identity', sequence='aa', cutoff=0)
ir.tl.define_clonotypes(adata_clean, receptor_arms='all', dual_ir='primary_only')

n_clonotypes = adata_clean.obs['clone_id'].nunique()
print(f'Defined {n_clonotypes} unique clonotypes')

# Clonal expansion analysis
# Categorizes cells by clone size
ir.tl.clonal_expansion(adata_clean)

# Expansion categories: 1 (singleton), 2, 3-10, >10
print('\nClonal expansion:')
print(adata_clean.obs['clonal_expansion'].value_counts())

# Calculate percentage expanded (clone size > 1)
expanded = (adata_clean.obs['clonal_expansion'] != '1').sum()
pct_expanded = expanded / adata_clean.n_obs * 100
print(f'Expanded cells: {pct_expanded:.1f}%')

# Plot clonal expansion by cell type (if annotated)
if 'cell_type' in adata_clean.obs.columns:
    ir.pl.clonal_expansion(adata_clean, groupby='cell_type')

# V gene usage
# Summarize which V genes are most common
ir.pl.vdj_usage(
    adata_clean,
    vdj_cols=['IR_VDJ_1_v_call'],  # TRB V gene
    full_names=False,
    max_segments=20
)

# Spectratype - CDR3 length distribution
# Normal repertoire shows Gaussian-like distribution
ir.pl.spectratype(adata_clean, chain='TRB', target_col='clone_id')

# Integration with gene expression
# Find genes associated with clonal expansion
adata_clean.obs['is_expanded'] = (adata_clean.obs['clonal_expansion'] != '1').astype(str)

sc.tl.rank_genes_groups(adata_clean, groupby='is_expanded', method='wilcoxon')
print('\nTop genes in expanded vs singleton cells:')
markers = sc.get.rank_genes_groups_df(adata_clean, group='True').head(10)
print(markers[['names', 'logfoldchanges', 'pvals_adj']])

# Export clonotype table
clonotype_df = adata_clean.obs[[
    'clone_id', 'clonal_expansion',
    'IR_VDJ_1_junction_aa', 'IR_VJ_1_junction_aa',
    'IR_VDJ_1_v_call', 'IR_VDJ_1_j_call'
]].drop_duplicates()
clonotype_df.to_csv('clonotypes_summary.csv', index=False)

print('\nAnalysis complete!')
print(f'Results saved to clonotypes_summary.csv')
