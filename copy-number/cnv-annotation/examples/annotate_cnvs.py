#!/usr/bin/env python3
'''Annotate CNV calls with genes'''

import pybedtools
import pandas as pd
import sys

def annotate_cnvs(cnv_file, genes_file, output_file):
    '''Annotate CNVs with overlapping genes'''
    cnvs = pybedtools.BedTool(cnv_file)
    genes = pybedtools.BedTool(genes_file)

    annotated = cnvs.intersect(genes, wa=True, wb=True)

    results = []
    for feature in annotated:
        results.append({
            'cnv_chrom': feature[0],
            'cnv_start': feature[1],
            'cnv_end': feature[2],
            'cnv_type': feature[3] if len(feature) > 3 else '.',
            'gene_chrom': feature[4] if len(feature) > 4 else '.',
            'gene_start': feature[5] if len(feature) > 5 else '.',
            'gene_end': feature[6] if len(feature) > 6 else '.',
            'gene_name': feature[7] if len(feature) > 7 else '.'
        })

    df = pd.DataFrame(results)

    # Summarize by CNV
    summary = df.groupby(['cnv_chrom', 'cnv_start', 'cnv_end', 'cnv_type']).agg({
        'gene_name': lambda x: ','.join(sorted(set(x)))
    }).reset_index()
    summary.columns = ['chrom', 'start', 'end', 'type', 'genes']

    summary.to_csv(output_file, sep='\t', index=False)
    print(f'Annotated {len(summary)} CNVs')
    print(f'Output: {output_file}')

if __name__ == '__main__':
    cnv_file = sys.argv[1] if len(sys.argv) > 1 else 'cnvs.bed'
    genes_file = sys.argv[2] if len(sys.argv) > 2 else 'genes.bed'
    output_file = sys.argv[3] if len(sys.argv) > 3 else 'cnvs_annotated.tsv'
    annotate_cnvs(cnv_file, genes_file, output_file)
