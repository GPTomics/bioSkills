#!/usr/bin/env python3
'''
Parse Infernal cmscan --fmt 2 tabular output and summarize ncRNA family assignments.
--fmt 2 inserts an idx column (0) and a clan column (5) versus --fmt 1, so every field
index shifts; a fmt-1 parser silently reads the wrong columns (score and E-value especially).
'''
# Reference: biopython 1.83+, infernal 1.1.4+, pandas 2.2+ | Verify API if version differs

import pandas as pd
from pathlib import Path
from Bio import SeqIO


def parse_cmscan_tblout(tblout_file):
    '''
    Parse Infernal cmscan --tblout --fmt 2 output.

    0-based fmt-2 cmscan indices: idx 0, family name 1, family accession 2, seq name 3,
    seq accession 4, clan 5, mdl type 6, mdl_from 7, mdl_to 8, seq_from 9, seq_to 10,
    strand 11, trunc 12, pass 13, gc 14, bias 15, score 16, E-value 17, inc 18, olp 19.
    (In cmsearch the family and sequence names are swapped: seq name is column 1.)
    '''
    rows = []
    with open(tblout_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.split()
            if len(fields) < 20:
                continue
            rows.append({
                'family_name': fields[1], 'family_accession': fields[2],
                'seq_name': fields[3], 'clan': fields[5],
                'mdl_from': int(fields[7]), 'mdl_to': int(fields[8]),
                'seq_from': int(fields[9]), 'seq_to': int(fields[10]),
                'strand': fields[11], 'trunc': fields[12],
                'gc': float(fields[14]), 'bias': float(fields[15]),
                'score': float(fields[16]), 'evalue': float(fields[17]),
                'inc': fields[18], 'olp': fields[19],
                'description': ' '.join(fields[26:]) if len(fields) > 26 else ''
            })
    return pd.DataFrame(rows)


def deoverlap_clans(df):
    '''Drop hits marked '=' (dominated by a higher-scoring clanmate), as Rfam does.'''
    return df[df['olp'] != '='].copy()


def filter_hits(df, evalue_threshold=1e-5, min_score=None):
    '''Filter for significant hits (use --cut_ga at search time for Rfam; E-value is DB-size-dependent).'''
    filtered = df[df['evalue'] <= evalue_threshold].copy()
    if min_score is not None:
        filtered = filtered[filtered['score'] >= min_score]
    return filtered.sort_values('score', ascending=False).reset_index(drop=True)


def summarize_families(df):
    '''Count hits per ncRNA family.'''
    return df.groupby(['family_name', 'family_accession']).agg(
        count=('seq_name', 'count'), best_score=('score', 'max'), best_evalue=('evalue', 'min')
    ).sort_values('count', ascending=False).reset_index()


def extract_hit_sequences(fasta_file, hits_df, output_file):
    '''Extract the sequence of each hit (a '-' strand hit has seq_from > seq_to; reverse-complement it).'''
    seqs = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    records = []
    for _, hit in hits_df.iterrows():
        if hit['seq_name'] not in seqs:
            continue
        start, end = sorted([hit['seq_from'], hit['seq_to']])
        subseq = seqs[hit['seq_name']][start - 1:end]
        if hit['strand'] == '-':
            subseq = subseq.reverse_complement()
        subseq.id = f'{hit["seq_name"]}_{start}_{end}_{hit["family_name"]}'
        subseq.description = f'family={hit["family_name"]} score={hit["score"]:.1f} E={hit["evalue"]:.1e}'
        records.append(subseq)
    SeqIO.write(records, output_file, 'fasta')
    print(f'Extracted {len(records)} hit sequences to {output_file}')
    return records


if __name__ == '__main__':
    tblout = 'rfam_results.tbl'
    query_fasta = 'query.fa'

    if not Path(tblout).exists():
        print(f'No results file at {tblout}; run rfam_search.sh first.')
        exit(0)

    print('=== Parsing Infernal --fmt 2 output ===')
    df = deoverlap_clans(parse_cmscan_tblout(tblout))
    print(f'Hits after clan deoverlapping: {len(df)}')

    significant = filter_hits(df, evalue_threshold=1e-5)
    print(f'Significant hits (E < 1e-5): {len(significant)}')

    if len(significant) > 0:
        cols = ['family_name', 'seq_name', 'seq_from', 'seq_to', 'strand', 'score', 'evalue']
        print('\n=== Top hits ===')
        print(significant[cols].head(20).to_string(index=False))
        print('\n=== Family summary ===')
        print(summarize_families(significant).to_string(index=False))
        if Path(query_fasta).exists():
            print('\n=== Extracting hit sequences ===')
            extract_hit_sequences(query_fasta, significant, 'rfam_hits.fa')
