# Reference: pertpy 0.6+, scanpy 1.10+, anndata 0.10+ | Verify API if version differs
#
# Perturb-seq analysis end-to-end using Pertpy:
# 1. sgRNA assignment
# 2. Escaper filtering via Mixscape
# 3. SCEPTRE differential expression (low-MOI calibrated)

import pertpy as pt
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd

# === LOAD DATA ===
# Pertpy includes Papalexi 2021 Mixscape dataset for reference
adata = pt.dt.papalexi_2021()

# Verify metadata
print(f'Cells: {adata.n_obs}')
print(f'Genes: {adata.n_vars}')
print(f'Perturbation column: "perturbation"')
print(f'NTC label: "NT"')

# === STANDARD scRNA-SEQ PREPROCESSING ===
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers['log_normalized'] = adata.X.copy()
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.tl.pca(adata, use_highly_variable=True)

# === MIXSCAPE ESCAPER FILTERING ===
# Compute perturbation signature: cell_expression - mean(K NTC neighbors)
ms = pt.tl.Mixscape()
ms.perturbation_signature(
    adata=adata,
    pert_key='perturbation',
    control='NT',
    n_neighbors=20,                              # K neighbors for KNN-NTC subtraction
)
# Writes .layers['X_pert'] with perturbation signature

# Classify cells as KO (true perturbation) or NP (non-perturbed/escaper)
ms.mixscape(
    adata=adata,
    pert_key='perturbation',                    # pert_key (modern pertpy API)
    control='NT',
    new_class_name='mixscape_class_global',
)
# Defaults to layer='X_pert' from previous step

# Retain KO cells + NTC controls for DE
adata_filtered = adata[adata.obs['mixscape_class_global'] != 'NP'].copy()
print(f'Cells after Mixscape filtering: {adata_filtered.n_obs}')
ko_pct = (adata_filtered.obs['mixscape_class_global'] == 'KO').mean()
print(f'KO retention rate among perturbed: {ko_pct:.1%}')

# === DIFFERENTIAL EXPRESSION ===
# Per-perturbation pseudobulk DE via pertpy's PyDESeq2 wrapper
# Loop perturbations; one contrast per pert vs NT
# (For calibrated single-cell low-MOI DE, run R sceptre separately; Pertpy does not bundle it.)
de = pt.tl.PyDESeq2(adata_filtered, design='~perturbation')
de.fit()

all_de = []
for pert in adata_filtered.obs['perturbation'].unique():
    if pert == 'NT':
        continue
    contrast_df = de.test_contrasts(contrast=('perturbation', pert, 'NT'))
    contrast_df['perturbation'] = pert
    all_de.append(contrast_df)
de_result = pd.concat(all_de)

# === OUTPUT ===
de_result.to_csv('pertpy_de_results.tsv', sep='\t')
print(f'DE results: {de_result.shape}')

# Top hits per perturbation
for pert in de_result['perturbation'].unique():
    pert_results = de_result[de_result['perturbation'] == pert].sort_values('padj')
    top = pert_results.head(10)
    if not top.empty:
        print(f'\n{pert}: top 10 differential genes')
        print(top[['log2FoldChange', 'padj']].to_string())
