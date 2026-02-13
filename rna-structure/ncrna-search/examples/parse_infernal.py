#!/usr/bin/env python3
'''
Parse Infernal cmscan/cmsearch output and summarize ncRNA family assignments.
'''
# Reference: biopython 1.83+, infernal 1.1+, pandas 2.2+ | Verify API if version differs

import pandas as pd
from pathlib import Path
from Bio import SeqIO


def parse_cmscan_tblout(tblout_file):
    '''Parse Infernal --tblout --fmt 2 tabular output.'''
    rows = []
    with open(tblout_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split()
            if len(fields) < 18:
                continue
            rows.append({
                'target_name': fields[0],
                'target_accession': fields[1],
                'query_name': fields[2],
                'query_accession': fields[3],
                'mdl_type': fields[4],
                'mdl_from': int(fields[5]),
                'mdl_to': int(fields[6]),
                'seq_from': int(fields[7]),
                'seq_to': int(fields[8]),
                'strand': fields[9],
                'trunc': fields[10],
                'pass': fields[11],
                'gc': float(fields[12]),
                'bias': float(fields[13]),
                'score': float(fields[14]),
                'evalue': float(fields[15]),
                'inc': fields[16],
                'description': ' '.join(fields[17:])
            })
    return pd.DataFrame(rows)


def filter_hits(df, evalue_threshold=1e-5, min_score=None):
    '''Filter for significant hits.'''
    filtered = df[df['evalue'] <= evalue_threshold].copy()
    if min_score is not None:
        filtered = filtered[filtered['score'] >= min_score]
    return filtered.sort_values('score', ascending=False).reset_index(drop=True)


def summarize_families(df):
    '''Count hits per ncRNA family.'''
    summary = df.groupby(['target_name', 'target_accession']).agg(
        count=('query_name', 'count'),
        best_score=('score', 'max'),
        best_evalue=('evalue', 'min'),
        description=('description', 'first')
    ).sort_values('count', ascending=False).reset_index()
    return summary


def extract_hit_sequences(fasta_file, hits_df, output_file):
    '''Extract sequences corresponding to Infernal hits.'''
    seqs = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    records = []

    for _, hit in hits_df.iterrows():
        query = hit['query_name']
        if query not in seqs:
            continue

        start = min(hit['seq_from'], hit['seq_to'])
        end = max(hit['seq_from'], hit['seq_to'])

        subseq = seqs[query][start-1:end]
        if hit['strand'] == '-':
            subseq = subseq.reverse_complement()

        subseq.id = f'{query}_{start}_{end}_{hit["target_name"]}'
        subseq.description = f'family={hit["target_name"]} score={hit["score"]:.1f} E={hit["evalue"]:.1e}'
        records.append(subseq)

    SeqIO.write(records, output_file, 'fasta')
    print(f'Extracted {len(records)} hit sequences to {output_file}')
    return records


if __name__ == '__main__':
    # Example: parse cmscan output
    tblout = 'rfam_results.tbl'
    query_fasta = 'query.fa'

    if not Path(tblout).exists():
        print(f'No results file found at {tblout}')
        print('Run rfam_search.sh first or provide a valid .tbl file')
        exit(0)

    print('=== Parsing Infernal output ===')
    df = parse_cmscan_tblout(tblout)
    print(f'Total hits: {len(df)}')

    # E-value < 1e-5: high-confidence family assignment
    significant = filter_hits(df, evalue_threshold=1e-5)
    print(f'Significant hits (E < 1e-5): {len(significant)}')

    if len(significant) > 0:
        print('\n=== Top hits ===')
        cols = ['target_name', 'query_name', 'seq_from', 'seq_to', 'strand', 'score', 'evalue']
        print(significant[cols].head(20).to_string(index=False))

        print('\n=== Family summary ===')
        families = summarize_families(significant)
        print(families.to_string(index=False))

        if Path(query_fasta).exists():
            print('\n=== Extracting hit sequences ===')
            extract_hit_sequences(query_fasta, significant, 'rfam_hits.fa')
