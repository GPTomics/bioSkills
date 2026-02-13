#!/usr/bin/env python3
'''
Parse and combine non-coding RNA annotations from Infernal and tRNAscan-SE.
'''
# Reference: pandas 2.2+ | Verify API if version differs

import pandas as pd
from collections import defaultdict


def parse_infernal_tbl(tbl_file):
    '''Parse Infernal cmscan --tblout --fmt 2 output.'''
    hits = []
    with open(tbl_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 18:
                continue
            strand = '+' if parts[9] == '+' else '-'
            start, end = (int(parts[7]), int(parts[8])) if strand == '+' else (int(parts[8]), int(parts[7]))
            hits.append({
                'rfam_acc': parts[1],
                'rfam_name': parts[2],
                'seqid': parts[3],
                'start': start,
                'end': end,
                'strand': strand,
                'evalue': float(parts[16]),
                'score': float(parts[15]),
            })
    df = pd.DataFrame(hits)
    if len(df) > 0:
        df['rna_type'] = df['rfam_name'].apply(classify_rfam)
    return df


def classify_rfam(name):
    '''Classify Rfam family into broad ncRNA categories.'''
    n = name.lower()
    if any(k in n for k in ['rrna', 'ssu', 'lsu', '5s_rrna', '5_8s']):
        return 'rRNA'
    if 'trna' in n:
        return 'tRNA'
    if any(k in n for k in ['snora', 'snord', 'snorna', 'haca_box', 'cd_box']):
        return 'snoRNA'
    if 'mir' in n:
        return 'miRNA'
    if 'riboswitch' in n or 'thermoregulator' in n:
        return 'riboswitch'
    if any(k in n for k in ['ires', 'leader', 'utr']):
        return 'cis-reg'
    if 'snrna' in n or n.startswith('u') and n[1:].isdigit():
        return 'snRNA'
    return 'other_ncRNA'


def parse_trnascan_gff(gff_file):
    '''Parse tRNAscan-SE GFF3 output.'''
    trnas = []
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'tRNA':
                continue
            attrs = {}
            for item in parts[8].split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    attrs[k] = v
            trnas.append({
                'seqid': parts[0],
                'start': int(parts[3]),
                'end': int(parts[4]),
                'strand': parts[6],
                'score': float(parts[5]) if parts[5] != '.' else 0,
                'isotype': attrs.get('isotype', 'unknown'),
                'anticodon': attrs.get('anticodon', 'unknown'),
                'rna_type': 'tRNA',
            })
    return pd.DataFrame(trnas)


def combine_and_summarize(infernal_tbl, trnascan_gff):
    '''Combine ncRNA annotations, preferring tRNAscan-SE for tRNAs.'''
    infernal_df = parse_infernal_tbl(infernal_tbl)
    trna_df = parse_trnascan_gff(trnascan_gff)

    infernal_no_trna = infernal_df[infernal_df['rna_type'] != 'tRNA'] if len(infernal_df) > 0 else infernal_df

    summary = defaultdict(int)
    if len(infernal_no_trna) > 0:
        for rna_type, count in infernal_no_trna['rna_type'].value_counts().items():
            summary[rna_type] = count
    summary['tRNA'] = len(trna_df)

    print('=== Combined ncRNA Summary ===')
    total = 0
    for rna_type in sorted(summary.keys()):
        print(f'  {rna_type}: {summary[rna_type]}')
        total += summary[rna_type]
    print(f'  Total ncRNAs: {total}')

    if len(trna_df) > 0:
        print(f'\n=== tRNA Isotype Distribution ===')
        isotype_counts = trna_df['isotype'].value_counts()
        for isotype, count in isotype_counts.items():
            print(f'  {isotype}: {count}')

    if len(infernal_no_trna) > 0:
        print(f'\n=== Top Rfam Families (non-tRNA) ===')
        top_families = infernal_no_trna['rfam_name'].value_counts().head(15)
        for family, count in top_families.items():
            print(f'  {family}: {count}')

    return infernal_no_trna, trna_df, summary


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 3:
        print('Usage: parse_ncrna.py <infernal_hits.tbl> <trnascan.gff3>')
        sys.exit(1)

    infernal_df, trna_df, summary = combine_and_summarize(sys.argv[1], sys.argv[2])

    if len(infernal_df) > 0:
        infernal_df.to_csv('infernal_parsed.tsv', sep='\t', index=False)
    if len(trna_df) > 0:
        trna_df.to_csv('trnascan_parsed.tsv', sep='\t', index=False)
    print('\nParsed results saved to infernal_parsed.tsv and trnascan_parsed.tsv')
