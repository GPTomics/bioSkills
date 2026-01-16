#!/usr/bin/env python3
'''Parse GTF files and extract gene information using gtfparse.'''

import gtfparse
import pandas as pd

# For demo, create a minimal GTF-like structure
# In practice, use: df = gtfparse.read_gtf('annotation.gtf')
demo_data = {
    'seqname': ['chr1', 'chr1', 'chr1', 'chr1', 'chr17', 'chr17'],
    'source': ['HAVANA'] * 6,
    'feature': ['gene', 'transcript', 'exon', 'exon', 'gene', 'transcript'],
    'start': [11869, 11869, 11869, 12613, 7668402, 7668402],
    'end': [14409, 14409, 12227, 12721, 7687550, 7687550],
    'score': ['.'] * 6,
    'strand': ['+', '+', '+', '+', '-', '-'],
    'frame': ['.'] * 6,
    'gene_id': ['ENSG00000223972', 'ENSG00000223972', 'ENSG00000223972',
                'ENSG00000223972', 'ENSG00000141510', 'ENSG00000141510'],
    'gene_name': ['DDX11L1', 'DDX11L1', 'DDX11L1', 'DDX11L1', 'TP53', 'TP53'],
    'gene_type': ['transcribed_unprocessed_pseudogene'] * 4 + ['protein_coding'] * 2,
    'transcript_id': ['', 'ENST00000456328', 'ENST00000456328', 'ENST00000456328',
                      '', 'ENST00000269305'],
}
df = pd.DataFrame(demo_data)

print('=== GTF DataFrame Structure ===')
print(f'Total features: {len(df)}')
print(f'Columns: {list(df.columns)}')

print('\n=== Feature Counts ===')
print(df['feature'].value_counts())

print('\n=== Gene Information ===')
genes = df[df['feature'] == 'gene']
for _, gene in genes.iterrows():
    print(f"{gene['gene_name']} ({gene['gene_id']})")
    print(f"  Location: {gene['seqname']}:{gene['start']}-{gene['end']} ({gene['strand']})")
    print(f"  Type: {gene['gene_type']}")

print('\n=== Extract Protein-Coding Genes ===')
coding = df[(df['feature'] == 'gene') & (df['gene_type'] == 'protein_coding')]
print(f'Found {len(coding)} protein-coding genes')
for _, gene in coding.iterrows():
    print(f"  {gene['gene_name']}")

print('\n=== Convert to BED Format ===')
bed = genes[['seqname', 'start', 'end', 'gene_name', 'gene_id', 'strand']].copy()
bed['start'] = bed['start'] - 1  # Convert to 0-based
bed.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand']
print(bed.to_string(index=False))

print('\n=== Get Exons for DDX11L1 ===')
ddx11l1_exons = df[(df['gene_name'] == 'DDX11L1') & (df['feature'] == 'exon')]
for _, exon in ddx11l1_exons.iterrows():
    print(f"  Exon: {exon['seqname']}:{exon['start']}-{exon['end']}")

# Save to files
genes_bed = genes[['seqname', 'start', 'end', 'gene_name']].copy()
genes_bed['start'] = genes_bed['start'] - 1
genes_bed['score'] = 0
genes_bed['strand'] = genes['strand']
genes_bed = genes_bed[['seqname', 'start', 'end', 'gene_name', 'score', 'strand']]
genes_bed.to_csv('genes_from_gtf.bed', sep='\t', header=False, index=False)
print('\nSaved genes_from_gtf.bed')
