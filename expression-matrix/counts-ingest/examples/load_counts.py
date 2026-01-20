#!/usr/bin/env python3
'''Load count matrices from various formats'''

import pandas as pd
import anndata as ad
import numpy as np

def load_featurecounts(path):
    '''Load featureCounts output'''
    df = pd.read_csv(path, sep='\t', comment='#', index_col=0)
    df = df.iloc[:, 5:]  # Drop annotation columns
    df.columns = [c.replace('.bam', '') for c in df.columns]
    return df

def load_salmon_quant(quant_dirs):
    '''Load Salmon quantification from multiple samples'''
    from collections import defaultdict
    counts = defaultdict(dict)

    for qdir in quant_dirs:
        sample = qdir.split('/')[-1]
        quant = pd.read_csv(f'{qdir}/quant.sf', sep='\t', index_col=0)
        for gene, row in quant.iterrows():
            counts[gene][sample] = row['NumReads']

    return pd.DataFrame(counts).T

def load_rsem(paths):
    '''Load RSEM gene-level output'''
    dfs = []
    for path in paths:
        sample = path.split('/')[-1].replace('.genes.results', '')
        df = pd.read_csv(path, sep='\t', index_col=0)
        df = df[['expected_count']].rename(columns={'expected_count': sample})
        dfs.append(df)
    return pd.concat(dfs, axis=1)

def counts_to_anndata(counts_df, metadata_df=None):
    '''Convert count matrix to AnnData'''
    adata = ad.AnnData(X=counts_df.T.values,
                       obs=pd.DataFrame(index=counts_df.columns),
                       var=pd.DataFrame(index=counts_df.index))
    if metadata_df is not None:
        adata.obs = metadata_df.loc[adata.obs_names]
    return adata

if __name__ == '__main__':
    # Example usage
    print('Load featureCounts: load_featurecounts("counts.txt")')
    print('Load Salmon: load_salmon_quant(["sample1_quant", "sample2_quant"])')
