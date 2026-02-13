#!/usr/bin/env python3
'''
Merge eggNOG-mapper and InterProScan functional annotations.
'''
# Reference: pandas 2.2+ | Verify API if version differs

import pandas as pd
from pathlib import Path


def parse_eggnog(annotations_file):
    '''Parse eggNOG-mapper v2.1.x annotations output.'''
    df = pd.read_csv(annotations_file, sep='\t', comment='#', header=None, skiprows=5)
    col_names = [
        'query', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs',
        'max_annot_lvl', 'COG_category', 'Description', 'Preferred_name',
        'GOs', 'EC', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module',
        'KEGG_Reaction', 'KEGG_rclass', 'BRITE', 'KEGG_TC', 'CAZy',
        'BiGG_Reaction', 'PFAMs'
    ]
    df.columns = col_names[:len(df.columns)]
    return df


def parse_interproscan_tsv(tsv_file):
    '''Parse InterProScan TSV output.'''
    col_names = [
        'protein_id', 'md5', 'length', 'analysis', 'signature_acc',
        'signature_desc', 'start', 'stop', 'score', 'status', 'date',
        'interpro_acc', 'interpro_desc', 'go_terms', 'pathways'
    ]
    df = pd.read_csv(tsv_file, sep='\t', header=None, names=col_names)
    return df


def combine_go_terms(*go_strings):
    '''Combine GO terms from multiple sources, deduplicating.'''
    terms = set()
    for go_str in go_strings:
        if pd.notna(go_str) and go_str != '-':
            terms.update(t.strip() for t in str(go_str).replace('|', ',').split(',') if t.strip().startswith('GO:'))
    return ','.join(sorted(terms)) if terms else '-'


def merge_results(eggnog_file, interpro_file, output_file):
    '''Merge eggNOG-mapper and InterProScan results into a single table.'''
    eggnog_df = parse_eggnog(eggnog_file)

    interpro_df = parse_interproscan_tsv(interpro_file)
    interpro_summary = interpro_df.groupby('protein_id').agg({
        'signature_acc': lambda x: ';'.join(sorted(x.dropna().unique())),
        'interpro_acc': lambda x: ';'.join(sorted(x.dropna().astype(str).unique())),
        'go_terms': lambda x: '|'.join(x.dropna().unique()),
        'analysis': lambda x: ';'.join(sorted(x.unique())),
    }).reset_index()
    interpro_summary.columns = ['query', 'ipr_signatures', 'ipr_ids', 'ipr_go', 'ipr_databases']

    merged = eggnog_df.merge(interpro_summary, on='query', how='outer')
    merged['combined_go'] = merged.apply(
        lambda row: combine_go_terms(row.get('GOs', ''), row.get('ipr_go', '')), axis=1
    )

    merged.to_csv(output_file, sep='\t', index=False)
    return merged


def annotation_coverage(merged_df):
    '''Report annotation coverage statistics.'''
    total = len(merged_df)
    stats = {
        'total_proteins': total,
        'with_go': (merged_df['combined_go'] != '-').sum(),
        'with_kegg': merged_df['KEGG_ko'].notna().sum() if 'KEGG_ko' in merged_df else 0,
        'with_pfam': merged_df['PFAMs'].notna().sum() if 'PFAMs' in merged_df else 0,
        'with_ec': merged_df['EC'].notna().sum() if 'EC' in merged_df else 0,
        'with_interpro': merged_df['ipr_ids'].notna().sum() if 'ipr_ids' in merged_df else 0,
        'with_description': (merged_df['Description'].notna() & (merged_df['Description'] != '-')).sum() if 'Description' in merged_df else 0,
    }

    # At least one functional annotation from any source
    has_any = (
        (merged_df['combined_go'] != '-') |
        merged_df.get('PFAMs', pd.Series(dtype=str)).notna() |
        merged_df.get('KEGG_ko', pd.Series(dtype=str)).notna() |
        merged_df.get('ipr_ids', pd.Series(dtype=str)).notna()
    ).sum()
    stats['with_any_annotation'] = has_any

    print('=== Annotation Coverage ===')
    for key, val in stats.items():
        if key == 'total_proteins':
            print(f'  {key}: {val}')
        else:
            print(f'  {key}: {val} ({val/total:.1%})')

    return stats


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 3:
        print('Usage: merge_annotations.py <eggnog_annotations> <interpro_tsv> [output.tsv]')
        sys.exit(1)

    eggnog_file = sys.argv[1]
    interpro_file = sys.argv[2]
    output_file = sys.argv[3] if len(sys.argv) > 3 else 'merged_annotations.tsv'

    merged = merge_results(eggnog_file, interpro_file, output_file)
    annotation_coverage(merged)
    print(f'\nMerged annotations saved to {output_file}')
