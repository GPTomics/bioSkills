#!/usr/bin/env python3
'''Annotate CNV segments with genes, overlap extent, and dosage sensitivity.

Goes beyond raw gene overlap: records the fraction of each gene covered (whole-gene
vs partial), joins ClinGen haploinsufficiency/triplosensitivity scores, and flags
direction-consistent driver hits (oncogene amplified, tumor suppressor deleted).
Overlap is not consequence -- this script produces the evidence to decide it.
'''
# Reference: bedtools 2.31+, pybedtools 0.9+, pandas 2.2+ | Verify API if version differs

import sys
import pybedtools
import pandas as pd

GAIN_LOG2 = 0.3   # visual/triage cut; not a calling threshold
LOSS_LOG2 = -0.3


def annotate_cnvs(cnv_bed, genes_bed, clingen_dosage=None, cgc_file=None,
                  output='cnv_annotated.tsv'):
    '''Annotate CNV segments. cnv_bed columns: chrom/start/end/log2.
    genes_bed columns: chrom/start/end/gene.'''
    cnvs = pybedtools.BedTool(cnv_bed)
    genes = pybedtools.BedTool(genes_bed)

    # -wo appends the overlap length so partial vs whole-gene overlap is recoverable.
    rows = []
    for f in cnvs.intersect(genes, wo=True):
        rows.append({
            'chrom': f[0], 'start': int(f[1]), 'end': int(f[2]), 'log2': float(f[3]),
            'gene': f[7], 'gene_start': int(f[5]), 'gene_end': int(f[6]),
            'overlap_bp': int(f[-1]),
        })
    df = pd.DataFrame(rows)
    if df.empty:
        print('No CNV-gene overlaps found; check genome builds match.')
        return df

    df['gene_len'] = df['gene_end'] - df['gene_start']
    df['gene_frac_covered'] = (df['overlap_bp'] / df['gene_len']).clip(upper=1.0)
    df['whole_gene'] = df['gene_frac_covered'] >= 0.99

    if clingen_dosage:
        dosage = pd.read_csv(clingen_dosage, sep='\t')  # gene, HI_score, TS_score
        df = df.merge(dosage, on='gene', how='left')

    if cgc_file:
        cgc = pd.read_csv(cgc_file, sep='\t')
        role = dict(zip(cgc['Gene Symbol'], cgc['Role in Cancer'].fillna('')))
        df['driver_role'] = df['gene'].map(role).fillna('')
        df['driver_hit'] = (
            ((df['log2'] > GAIN_LOG2) & df['driver_role'].str.contains('oncogene')) |
            ((df['log2'] < LOSS_LOG2) & df['driver_role'].str.contains('TSG')))

    df.to_csv(output, sep='\t', index=False)
    print(f'Annotated {df["gene"].nunique()} genes across '
          f'{df[["chrom", "start", "end"]].drop_duplicates().shape[0]} CNV segments')
    print(f'Output: {output}')
    return df


if __name__ == '__main__':
    cnv = sys.argv[1] if len(sys.argv) > 1 else 'cnvs.bed'
    genes = sys.argv[2] if len(sys.argv) > 2 else 'genes.bed'
    dosage = sys.argv[3] if len(sys.argv) > 3 else None
    cgc = sys.argv[4] if len(sys.argv) > 4 else None
    annotate_cnvs(cnv, genes, dosage, cgc)
