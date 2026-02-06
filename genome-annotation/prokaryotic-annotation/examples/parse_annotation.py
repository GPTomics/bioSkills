#!/usr/bin/env python3
'''
Parse prokaryotic genome annotation GFF3 output and compute QC metrics.
'''

import gffutils
import pandas as pd
from Bio import SeqIO
from pathlib import Path


def load_gff_database(gff_file):
    '''Load GFF3 file into a queryable gffutils database.'''
    return gffutils.create_db(str(gff_file), ':memory:', merge_strategy='merge', sort_attribute_values=True)


def extract_cds_info(db):
    '''Extract CDS features with annotations into a DataFrame.'''
    records = []
    for cds in db.features_of_type('CDS'):
        records.append({
            'locus_tag': cds.attributes.get('locus_tag', [''])[0],
            'product': cds.attributes.get('product', ['unknown'])[0],
            'seqid': cds.seqid,
            'start': cds.start,
            'end': cds.end,
            'strand': cds.strand,
            'length_bp': cds.end - cds.start + 1
        })
    return pd.DataFrame(records)


def compute_coding_density(db, genome_fasta):
    '''Compute coding density from GFF3 and genome FASTA.

    Typical prokaryotic coding density: 85-95%.
    Below 80% may indicate annotation gaps or many pseudogenes.
    Above 95% may indicate overlapping annotations needing review.
    '''
    genome_length = sum(len(rec.seq) for rec in SeqIO.parse(genome_fasta, 'fasta'))
    coding_bp = sum(cds.end - cds.start + 1 for cds in db.features_of_type('CDS'))
    return coding_bp / genome_length, coding_bp, genome_length


def feature_summary(db):
    '''Summarize all feature types in the annotation.'''
    feature_counts = {}
    for feature in db.all_features():
        ftype = feature.featuretype
        feature_counts[ftype] = feature_counts.get(ftype, 0) + 1
    return dict(sorted(feature_counts.items(), key=lambda x: -x[1]))


def annotation_qc_report(gff_file, genome_fasta):
    '''Generate a QC report for a prokaryotic annotation.'''
    db = load_gff_database(gff_file)

    print('=== Feature Summary ===')
    for ftype, count in feature_summary(db).items():
        print(f'  {ftype}: {count}')

    cds_df = extract_cds_info(db)
    print(f'\n=== CDS Statistics ===')
    print(f'  Total CDSs: {len(cds_df)}')
    print(f'  Median CDS length: {cds_df["length_bp"].median():.0f} bp')
    print(f'  Mean CDS length: {cds_df["length_bp"].mean():.0f} bp')

    hypothetical = cds_df[cds_df['product'].str.contains('hypothetical', case=False)]
    hypo_frac = len(hypothetical) / len(cds_df) if len(cds_df) > 0 else 0
    print(f'  Hypothetical proteins: {len(hypothetical)} ({hypo_frac:.1%})')

    density, coding_bp, genome_bp = compute_coding_density(db, genome_fasta)
    print(f'\n=== Coding Density ===')
    print(f'  Genome length: {genome_bp:,} bp')
    print(f'  Coding bases: {coding_bp:,} bp')
    print(f'  Coding density: {density:.1%}')

    if density < 0.80:
        print('  WARNING: Low coding density (<80%). Check assembly completeness.')
    elif density > 0.95:
        print('  WARNING: High coding density (>95%). Check for overlapping annotations.')

    return cds_df


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 3:
        print('Usage: parse_annotation.py <annotation.gff3> <genome.fasta>')
        sys.exit(1)

    cds_df = annotation_qc_report(sys.argv[1], sys.argv[2])
    cds_df.to_csv('cds_features.tsv', sep='\t', index=False)
    print('\nCDS features saved to cds_features.tsv')
