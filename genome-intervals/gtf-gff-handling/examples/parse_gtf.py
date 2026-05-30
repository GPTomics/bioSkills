#!/usr/bin/env python3
'''Walk the gene-model tree with gffutils and convert features to BED with the correct coordinate shift.'''
# Reference: gffutils 0.13+ | Verify API if version differs

import gffutils

GTF = '''chr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972"; gene_name "DDX11L1"; gene_biotype "lncRNA";
chr1\tHAVANA\ttranscript\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972"; transcript_id "ENST00000456328"; gene_name "DDX11L1";
chr1\tHAVANA\texon\t11869\t12227\t.\t+\t.\tgene_id "ENSG00000223972"; transcript_id "ENST00000456328"; exon_number "1";
chr1\tHAVANA\texon\t12613\t12721\t.\t+\t.\tgene_id "ENSG00000223972"; transcript_id "ENST00000456328"; exon_number "2";
chr1\tHAVANA\texon\t13221\t14409\t.\t+\t.\tgene_id "ENSG00000223972"; transcript_id "ENST00000456328"; exon_number "3";
'''

with open('demo.gtf', 'w') as fh:
    fh.write(GTF)

# This demo GTF already has gene+transcript lines, so disable inference (~100x faster on real GENCODE/Ensembl files)
db = gffutils.create_db('demo.gtf', ':memory:', force=True,
                        disable_infer_genes=True, disable_infer_transcripts=True,
                        merge_strategy='create_unique')

print('=== Feature counts ===')
for ft in db.featuretypes():
    print(ft, db.count_features_of_type(ft))

print('\n=== Walk gene -> transcript -> exon (gffutils keeps 1-based coords) ===')
gene = db['ENSG00000223972']
print(f'{gene.id} {gene.seqid}:{gene.start}-{gene.end} ({gene.strand})')
for tx in db.children(gene, featuretype=['mRNA', 'transcript'], order_by='start'):
    exons = list(db.children(tx, featuretype='exon', order_by='start'))
    introns = list(db.interfeatures(exons, new_featuretype='intron'))
    print(f'  {tx.id}: {len(exons)} exons, {len(introns)} derived introns')
    for intron in introns:
        print(f'    intron {intron.start}-{intron.end}')

print('\n=== Genes -> BED (subtract 1 from START only; end unchanged) ===')
bed_rows = [(g.seqid, g.start - 1, g.end, g.id, 0, g.strand) for g in db.features_of_type('gene')]
with open('genes_from_gtf.bed', 'w') as fh:
    for row in bed_rows:
        fh.write('\t'.join(str(x) for x in row) + '\n')
for row in bed_rows:
    print(row)
print('Saved genes_from_gtf.bed')
