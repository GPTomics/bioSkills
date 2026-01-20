#!/usr/bin/env python3
'''Map between gene ID systems'''

import mygene
import pandas as pd

def map_ensembl_to_symbol(ensembl_ids):
    '''Map Ensembl IDs to gene symbols'''
    mg = mygene.MyGeneInfo()

    clean_ids = [g.split('.')[0] for g in ensembl_ids]

    results = mg.querymany(clean_ids, scopes='ensembl.gene',
                           fields='symbol', species='human')

    mapping = {}
    for r in results:
        if 'symbol' in r:
            mapping[r['query']] = r['symbol']

    return mapping

def map_symbol_to_entrez(symbols):
    '''Map gene symbols to Entrez IDs'''
    mg = mygene.MyGeneInfo()

    results = mg.querymany(symbols, scopes='symbol',
                           fields='entrezgene', species='human')

    mapping = {}
    for r in results:
        if 'entrezgene' in r:
            mapping[r['query']] = r['entrezgene']

    return mapping

def add_gene_symbols(counts_df, id_column='index'):
    '''Add gene symbols to count matrix'''
    if id_column == 'index':
        ids = counts_df.index.tolist()
    else:
        ids = counts_df[id_column].tolist()

    mapping = map_ensembl_to_symbol(ids)

    if id_column == 'index':
        counts_df['symbol'] = counts_df.index.map(lambda x: mapping.get(x.split('.')[0], x))
    else:
        counts_df['symbol'] = counts_df[id_column].map(lambda x: mapping.get(x.split('.')[0], x))

    return counts_df

if __name__ == '__main__':
    test_ids = ['ENSG00000139618', 'ENSG00000141510', 'ENSG00000171862']
    mapping = map_ensembl_to_symbol(test_ids)
    print('Ensembl to Symbol mapping:')
    for k, v in mapping.items():
        print(f'  {k} -> {v}')
