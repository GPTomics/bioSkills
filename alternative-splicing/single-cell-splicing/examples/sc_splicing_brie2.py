#!/usr/bin/env python3
'''
Single-cell splicing analysis using BRIE2.
Estimates PSI values with uncertainty quantification for sparse scRNA-seq data.
'''

import subprocess
from pathlib import Path


def prepare_splicing_events(gtf_file, output_file):
    '''
    Generate splicing event annotations from GTF for BRIE2.
    '''
    import brie

    brie.preprocessing.get_events(
        gtf_file=gtf_file,
        out_file=output_file
    )

    print(f'Splicing events written to: {output_file}')


def count_splicing_reads(bam_file, gff_file, output_dir, barcode_file=None):
    '''
    Count reads for splicing events from 10X BAM.

    Args:
        bam_file: Possorted BAM from Cell Ranger
        gff_file: Splicing events GFF from prepare_splicing_events
        output_dir: Output directory for count matrices
        barcode_file: Optional filtered barcodes (barcodes.tsv)
    '''
    import brie

    kwargs = {
        'bam_file': bam_file,
        'gff_file': gff_file,
        'out_dir': output_dir
    }

    if barcode_file:
        kwargs['cell_file'] = barcode_file

    brie.preprocessing.count(**kwargs)

    print(f'BRIE counts written to: {output_dir}')


def run_brie2_inference(count_dir, n_epochs=400, batch_size=512):
    '''
    Run BRIE2 variational inference for PSI estimation.

    Returns AnnData with PSI estimates in layers.
    '''
    import brie

    adata = brie.read_h5ad(f'{count_dir}/brie_count.h5ad')

    # Run BRIE2 model
    # Probabilistic PSI estimation handles sparse data
    brie.fit(
        adata,
        layer='raw',
        n_epochs=n_epochs,
        batch_size=batch_size
    )

    # PSI estimates: adata.layers['Psi']
    # Uncertainty: adata.layers['Psi_var']

    return adata


def find_variable_splicing(adata_splice, cell_type_col='cell_type', min_cells=50):
    '''
    Identify splicing events variable across cell types.

    Args:
        adata_splice: AnnData from BRIE2 with PSI in layers
        cell_type_col: Column in obs with cell type labels
        min_cells: Minimum cells per type for reliable estimates
    '''
    import numpy as np
    import pandas as pd

    psi = adata_splice.layers['Psi']
    cell_types = adata_splice.obs[cell_type_col].unique()

    mean_psi = pd.DataFrame(index=adata_splice.var_names)

    for ct in cell_types:
        mask = adata_splice.obs[cell_type_col] == ct
        if mask.sum() >= min_cells:
            mean_psi[ct] = np.nanmean(psi[mask, :], axis=0)

    # Calculate variance across cell types
    psi_variance = mean_psi.var(axis=1)
    variable_events = psi_variance.sort_values(ascending=False)

    return variable_events, mean_psi


def pseudobulk_by_celltype(adata, cell_type_col='cell_type'):
    '''
    Create pseudobulk samples by aggregating cells within each type.
    Useful for running bulk splicing tools on scRNA-seq.
    '''
    import numpy as np
    import pandas as pd

    pseudobulk = {}

    for ct in adata.obs[cell_type_col].unique():
        mask = adata.obs[cell_type_col] == ct
        if mask.sum() >= 20:
            pseudobulk[ct] = np.array(adata.X[mask].sum(axis=0)).flatten()

    return pd.DataFrame(pseudobulk, index=adata.var_names)


def differential_splicing_pseudobulk(adata_splice, group1_cells, group2_cells):
    '''
    Compare splicing between two cell populations using pseudobulk approach.

    Args:
        adata_splice: AnnData with PSI estimates
        group1_cells: Boolean mask or list of cell barcodes for group 1
        group2_cells: Boolean mask or list of cell barcodes for group 2
    '''
    import numpy as np
    from scipy import stats

    psi = adata_splice.layers['Psi']

    if isinstance(group1_cells, list):
        mask1 = adata_splice.obs_names.isin(group1_cells)
        mask2 = adata_splice.obs_names.isin(group2_cells)
    else:
        mask1 = group1_cells
        mask2 = group2_cells

    results = []

    for i, event in enumerate(adata_splice.var_names):
        psi1 = psi[mask1, i]
        psi2 = psi[mask2, i]

        # Remove NaN values
        psi1 = psi1[~np.isnan(psi1)]
        psi2 = psi2[~np.isnan(psi2)]

        if len(psi1) >= 20 and len(psi2) >= 20:
            # Mann-Whitney U test
            stat, pval = stats.mannwhitneyu(psi1, psi2, alternative='two-sided')
            delta_psi = np.mean(psi1) - np.mean(psi2)

            results.append({
                'event': event,
                'mean_psi_group1': np.mean(psi1),
                'mean_psi_group2': np.mean(psi2),
                'delta_psi': delta_psi,
                'pvalue': pval,
                'n_cells_group1': len(psi1),
                'n_cells_group2': len(psi2)
            })

    import pandas as pd
    results_df = pd.DataFrame(results)

    # FDR correction
    from statsmodels.stats.multitest import multipletests
    if len(results_df) > 0:
        _, results_df['fdr'], _, _ = multipletests(results_df['pvalue'], method='fdr_bh')
        results_df = results_df.sort_values('fdr')

    return results_df


if __name__ == '__main__':
    # Example workflow
    print('BRIE2 Single-Cell Splicing Analysis')
    print('=' * 40)

    # Step 1: Prepare events
    # prepare_splicing_events('annotation.gtf', 'splicing_events.gff3')

    # Step 2: Count reads
    # count_splicing_reads(
    #     'possorted_genome_bam.bam',
    #     'splicing_events.gff3',
    #     'brie_counts/',
    #     'filtered_barcodes.tsv'
    # )

    # Step 3: Run BRIE2 inference
    # adata_splice = run_brie2_inference('brie_counts/')

    # Step 4: Find variable events
    # variable_events, mean_psi = find_variable_splicing(adata_splice)
    # print('Top variable splicing events:')
    # print(variable_events.head(20))

    print('Provide input files to run single-cell splicing analysis')
