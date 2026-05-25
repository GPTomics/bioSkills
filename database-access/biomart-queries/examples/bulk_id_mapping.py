'''Bulk ID mapping via pybiomart: Ensembl Gene -> HGNC + RefSeq + UniProt in one query.'''
# Reference: pybiomart 0.9+, Ensembl release 110+ | Verify API if version differs
from pybiomart import Server
import pandas as pd

server = Server(host='http://www.ensembl.org')
mart = server['ENSEMBL_MART_ENSEMBL']
ds = mart['hsapiens_gene_ensembl']


ensembl_ids = [
    'ENSG00000139618',  # BRCA2
    'ENSG00000141510',  # TP53
    'ENSG00000171862',  # PTEN
    'ENSG00000146648',  # EGFR
    'ENSG00000136997',  # MYC
]

print('=== Bulk ID mapping: Ensembl Gene -> HGNC + RefSeq + UniProt ===')
df = ds.query(
    attributes=[
        'ensembl_gene_id',
        'external_gene_name',
        'hgnc_id',
        'entrezgene_id',
        'refseq_mrna',
        'uniprotswissprot',
    ],
    filters={'ensembl_gene_id': ensembl_ids},
)
print(f'  Rows: {len(df)} (note: many-to-many cross-ref joins multiply rows)')
print(df.head(10).to_string(index=False))

print('\n=== Collapse to one row per gene (most common downstream pattern) ===')
collapsed = (df.groupby('Gene stable ID')
               .agg({'Gene name': 'first',
                     'HGNC ID': 'first',
                     'NCBI gene (formerly Entrezgene) ID': 'first',
                     'RefSeq mRNA ID': lambda x: ';'.join(set(filter(None, x))),
                     'UniProtKB/Swiss-Prot ID': lambda x: ';'.join(set(filter(None, x))) })
               .reset_index())
print(collapsed.to_string(index=False))


print('\n=== Discover attribute names if column labels differ ===')
print('First 10 of', len(ds.attributes), 'attributes:')
for a, info in list(ds.attributes.items())[:10]:
    print(f'  {a:<35} {info.display_name}')
