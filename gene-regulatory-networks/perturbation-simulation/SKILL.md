---
name: bio-gene-regulatory-networks-perturbation-simulation
description: Simulate transcription factor perturbation effects on cell state using CellOracle, which integrates GRN inference with in silico knockout and overexpression modeling. Predicts cell identity shifts and differentiation trajectory changes from TF perturbations. Use when predicting the effect of transcription factor knockouts, planning perturbation experiments, or identifying driver TFs for cell fate transitions.
tool_type: python
primary_tool: CellOracle
---

# Perturbation Simulation

Simulate transcription factor perturbation effects on cell state using CellOracle. Integrates GRN inference from scRNA-seq with base GRN from chromatin accessibility to predict cell identity shifts from TF knockouts or overexpression.

## CellOracle Overview

CellOracle constructs a GRN by combining:
1. A **base GRN** from accessible chromatin regions + motif scanning (defines possible TF-target links)
2. **scRNA-seq expression** data (learns active regulatory weights)

The base GRN can come from scATAC-seq, bulk ATAC-seq, or published chromatin data. CellOracle does NOT require paired multiome data -- any source of accessible regions works.

## Installation

```bash
pip install celloracle
```

## Step 1: Base GRN from Accessible Regions

### From scATAC-seq Peaks

```python
import celloracle as co
import pandas as pd
import numpy as np

# Load peak data (BED format: chr, start, end)
peaks = pd.read_csv('atac_peaks.bed', sep='\t', header=None, names=['chr', 'start', 'end'])

# Scan peaks for TF binding motifs using CellOracle's built-in scanner
# Uses gimmemotifs internally
tfi = co.motif_analysis.TFinfo(peak_data_frame=peaks, ref_genome='hg38')

# Scan for motifs
tfi.scan(fpr=0.02)  # false positive rate for motif matching

# Filter and format as base GRN
tfi.filter_motifs_by_score(threshold=10)
tfi.make_TFinfo_dataframe_and_target_gene_dataframe()

base_grn = tfi.to_dataframe()
base_grn.to_parquet('base_grn.parquet')
print(f'Base GRN: {len(base_grn)} TF-target links')
```

### From Published Chromatin Data

```python
# CellOracle provides pre-built base GRNs for common cell types
# Download from: https://github.com/morris-lab/CellOracle/wiki/
base_grn = co.data.load_mouse_scATAC_atlas_base_GRN(
    organism='Mouse',
    tissue='whole_brain'
)
```

## Step 2: GRN Construction from scRNA-seq

```python
import scanpy as sc
import celloracle as co

adata = sc.read_h5ad('clustered.h5ad')

# Ensure data is preprocessed: normalized, log-transformed, with PCA and clustering
oracle = co.Oracle()
oracle.import_anndata_as_raw_count(
    adata=adata,
    cluster_column_name='cell_type',
    embedding_name='X_umap'
)

# Load base GRN
base_grn = pd.read_parquet('base_grn.parquet')
oracle.import_TF_data(TF_info_matrix=base_grn)

# Fit GRN models per cluster
# Uses regularized linear regression (Bayesian Ridge) to learn TF-target weights
oracle.perform_PCA()
oracle.knn_imputation(n_pnn=30, balanced=True, b_sight=3000, b_maxl=1500)

# Fit GRN for all clusters
links = oracle.get_links(cluster_name_for_GRN_unit='cell_type', alpha=10, verbose_level=0)

# Filter links by statistical significance
# p-value threshold for keeping TF-target connections
links.filter_links(p=0.001, weight='coef_abs', threshold_number=2000)

# Inspect top regulatory connections
links.links_dict['T_cell'].sort_values('coef_abs', ascending=False).head(20)
```

## Step 3: Perturbation Simulation

### Knockout Simulation

```python
# Simulate TF knockout (set expression to 0)
oracle.simulate_shift(perturb_condition={'GATA1': 0.0}, n_propagation=3)

# Get perturbation scores
oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
oracle.calculate_embedding_shift(sigma_corr=0.05)

# Perturbation score: magnitude of predicted cell state shift
# Higher score = more affected by the perturbation
perturbation_scores = oracle.adata.obsm['delta_embedding']
shift_magnitude = np.sqrt((perturbation_scores ** 2).sum(axis=1))
oracle.adata.obs['GATA1_KO_shift'] = shift_magnitude
```

### Overexpression Simulation

```python
# Simulate TF overexpression (set to high value)
# Value is relative to max observed expression
oracle.simulate_shift(perturb_condition={'PAX5': 3.0}, n_propagation=3)
oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
oracle.calculate_embedding_shift(sigma_corr=0.05)
```

### Multi-TF Perturbation

```python
# Simulate multiple TF perturbations simultaneously
oracle.simulate_shift(
    perturb_condition={'GATA1': 0.0, 'SPI1': 0.0},  # double knockout
    n_propagation=3
)
oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
oracle.calculate_embedding_shift(sigma_corr=0.05)
```

## Visualization

### Quiver Plot (Vector Field)

```python
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(8, 8))

# Quiver plot shows predicted direction of cell state change
oracle.plot_quiver(
    ax=ax, scale=30,
    color=oracle.adata.obs['cell_type'],
    plot_whole_cells=True
)
ax.set_title('GATA1 KO - predicted cell state shifts')
plt.savefig('gata1_ko_quiver.pdf', bbox_inches='tight')
```

### Gradient Plot

```python
fig, ax = plt.subplots(1, 2, figsize=(16, 8))

# Perturbation score on embedding
sc.pl.embedding(oracle.adata, basis='umap', color='GATA1_KO_shift',
                cmap='Reds', ax=ax[0], show=False, title='Shift magnitude')

# Cell type reference
sc.pl.embedding(oracle.adata, basis='umap', color='cell_type',
                ax=ax[1], show=False, title='Cell types')

plt.savefig('gata1_ko_gradient.pdf', bbox_inches='tight')
```

### Systematic TF Screen

```python
# Screen multiple TFs to find drivers of cell fate
tfs_to_screen = ['GATA1', 'SPI1', 'CEBPA', 'PAX5', 'TCF7', 'RUNX1']

results = {}
for tf in tfs_to_screen:
    oracle.simulate_shift(perturb_condition={tf: 0.0}, n_propagation=3)
    oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
    oracle.calculate_embedding_shift(sigma_corr=0.05)

    shift = np.sqrt((oracle.adata.obsm['delta_embedding'] ** 2).sum(axis=1))
    results[tf] = {
        'mean_shift': shift.mean(),
        'max_shift': shift.max(),
        'affected_cells': (shift > shift.quantile(0.9)).sum()
    }

screen_df = pd.DataFrame(results).T.sort_values('mean_shift', ascending=False)
print(screen_df)
```

## Parameter Reference

| Parameter | Default | Description |
|-----------|---------|-------------|
| n_propagation | 3 | Signal propagation steps in GRN; higher = longer-range effects |
| n_neighbors | 200 | Neighbors for transition probability; adjust with dataset size |
| sigma_corr | 0.05 | Smoothing for embedding shift; lower = sharper gradients |
| alpha (GRN fit) | 10 | Regularization strength; higher = sparser GRN |
| p (link filter) | 0.001 | P-value cutoff for significant TF-target links |

## Related Skills

- scenic-regulons - TF regulon inference from scRNA-seq with pySCENIC
- multiomics-grn - Enhancer-driven GRNs with SCENIC+
- single-cell/trajectory-inference - Trajectory analysis for cell fate context
- single-cell/perturb-seq - Experimental perturbation data analysis
