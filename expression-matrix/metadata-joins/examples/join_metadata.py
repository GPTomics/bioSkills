#!/usr/bin/env python3
'''Join sample metadata with expression data'''

import pandas as pd
import anndata as ad

def join_metadata(counts_df, metadata_df, sample_col='sample_id'):
    '''Join metadata to count matrix columns'''
    if sample_col in metadata_df.columns:
        metadata_df = metadata_df.set_index(sample_col)

    missing = set(counts_df.columns) - set(metadata_df.index)
    if missing:
        print(f'Warning: {len(missing)} samples not in metadata')

    common = list(set(counts_df.columns) & set(metadata_df.index))
    counts_filtered = counts_df[common]
    metadata_filtered = metadata_df.loc[common]

    return counts_filtered, metadata_filtered

def create_deseq_input(counts_df, metadata_df, output_prefix='deseq'):
    '''Create files ready for DESeq2'''
    counts_df.to_csv(f'{output_prefix}_counts.csv')
    metadata_df.to_csv(f'{output_prefix}_metadata.csv')
    print(f'Saved: {output_prefix}_counts.csv, {output_prefix}_metadata.csv')

def add_metadata_to_anndata(adata, metadata_df, sample_col='sample_id'):
    '''Add metadata to AnnData object'''
    if sample_col in metadata_df.columns:
        metadata_df = metadata_df.set_index(sample_col)

    for col in metadata_df.columns:
        adata.obs[col] = metadata_df.loc[adata.obs_names, col]

    return adata

if __name__ == '__main__':
    # Example
    counts = pd.DataFrame({'sample1': [10, 20], 'sample2': [15, 25]}, index=['geneA', 'geneB'])
    meta = pd.DataFrame({'sample_id': ['sample1', 'sample2'], 'condition': ['ctrl', 'treat']})

    counts_j, meta_j = join_metadata(counts, meta)
    print('Joined successfully!')
