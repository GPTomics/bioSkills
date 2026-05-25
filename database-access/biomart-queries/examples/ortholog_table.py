'''Build a multi-species ortholog wide-table from one BioMart query; filter to 1:1 orthologs.'''
# Reference: pybiomart 0.9+, Ensembl release 110+ | Verify API if version differs
from pybiomart import Server
import pandas as pd

server = Server(host='http://www.ensembl.org')
mart = server['ENSEMBL_MART_ENSEMBL']
ds = mart['hsapiens_gene_ensembl']


print('=== Ortholog wide-table: human + mouse + zebrafish (chr17 subset) ===')
df = ds.query(
    attributes=[
        'ensembl_gene_id',
        'external_gene_name',
        'mmusculus_homolog_ensembl_gene',
        'mmusculus_homolog_orthology_type',
        'drerio_homolog_ensembl_gene',
        'drerio_homolog_orthology_type',
    ],
    filters={'chromosome_name': '17'},
)
print(f'  All chr17 genes with any ortholog row: {len(df)}')
print(df.head(8).to_string(index=False))


print('\n=== 1:1 across all three species ===')
# Column labels are the BioMart display names; check df.columns to confirm.
mouse_col = next(c for c in df.columns if 'Mouse' in c and 'type' in c)
zebra_col = next(c for c in df.columns if 'Zebrafish' in c and 'type' in c)
one2one = df[(df[mouse_col] == 'ortholog_one2one') &
             (df[zebra_col] == 'ortholog_one2one')]
print(f'  1:1 in mouse AND zebrafish: {len(one2one)}')
print(one2one.head(10).to_string(index=False))


print('\n=== Compare to per-gene Ensembl REST cost ===')
print(f'  {len(df)} chr17 genes via BioMart: one query, no rate limit')
print(f'  Same via Ensembl REST /homology/symbol: {len(df)} calls * 0.07s sleep = '
      f'{len(df) * 0.07 / 60:.1f} min minimum + HTTP overhead')
print(f'  Use REST only for <100 genes or real-time queries.')
