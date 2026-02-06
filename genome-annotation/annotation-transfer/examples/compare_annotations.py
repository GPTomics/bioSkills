#!/usr/bin/env python3
'''
Compare original and transferred annotations for quality assessment.
'''

import gffutils
import pandas as pd
from Bio import SeqIO


def load_gff_db(gff_file):
    return gffutils.create_db(str(gff_file), ':memory:', merge_strategy='merge', sort_attribute_values=True)


def annotation_stats(db, label='annotation'):
    '''Compute basic annotation statistics.'''
    genes = list(db.features_of_type('gene'))
    mrnas = list(db.features_of_type(['mRNA', 'transcript']))

    exon_counts = []
    gene_lengths = []
    for gene in genes:
        gene_lengths.append(gene.end - gene.start + 1)
        for tx in db.children(gene, featuretype=['mRNA', 'transcript']):
            exons = list(db.children(tx, featuretype='exon'))
            exon_counts.append(len(exons))

    stats = {
        'genes': len(genes),
        'transcripts': len(mrnas),
        'median_gene_length': pd.Series(gene_lengths).median() if gene_lengths else 0,
        'median_exons': pd.Series(exon_counts).median() if exon_counts else 0,
        'single_exon_pct': sum(1 for e in exon_counts if e == 1) / len(exon_counts) * 100 if exon_counts else 0,
    }

    print(f'\n=== {label} ===')
    for key, val in stats.items():
        if isinstance(val, float):
            print(f'  {key}: {val:.1f}')
        else:
            print(f'  {key}: {val}')

    return stats


def compare_annotations(reference_gff, transferred_gff):
    '''Compare reference and transferred annotations.'''
    ref_db = load_gff_db(reference_gff)
    tgt_db = load_gff_db(transferred_gff)

    ref_stats = annotation_stats(ref_db, 'Reference')
    tgt_stats = annotation_stats(tgt_db, 'Transferred')

    transfer_rate = tgt_stats['genes'] / ref_stats['genes'] if ref_stats['genes'] > 0 else 0

    print(f'\n=== Comparison ===')
    print(f'Gene transfer rate: {transfer_rate:.1%}')
    print(f'Transcript transfer rate: {tgt_stats["transcripts"] / ref_stats["transcripts"]:.1%}' if ref_stats['transcripts'] > 0 else '')

    # Transfer rate thresholds:
    # >95%: Excellent (same species, high-quality assemblies)
    # >80%: Good (same species, moderate assembly quality)
    # >60%: Acceptable (closely related species)
    # <60%: Poor (distant species or assembly issues)
    if transfer_rate > 0.95:
        print('Assessment: Excellent transfer quality')
    elif transfer_rate > 0.80:
        print('Assessment: Good transfer quality')
    elif transfer_rate > 0.60:
        print('Assessment: Acceptable, consider complementing with de novo prediction')
    else:
        print('Assessment: Low transfer rate, de novo prediction recommended')

    return ref_stats, tgt_stats


def validate_orfs(transferred_gff, target_fasta):
    '''Check ORF integrity of transferred gene models.'''
    genome = SeqIO.to_dict(SeqIO.parse(target_fasta, 'fasta'))
    db = load_gff_db(transferred_gff)

    results = {'valid': 0, 'no_start': 0, 'no_stop': 0, 'internal_stop': 0, 'frameshift': 0}
    total = 0

    for mrna in db.features_of_type(['mRNA', 'transcript']):
        cds_features = sorted(db.children(mrna, featuretype='CDS'), key=lambda x: x.start)
        if not cds_features:
            continue

        total += 1
        cds_seq = ''
        for cds in cds_features:
            seq = str(genome[cds.seqid].seq[cds.start - 1:cds.end])
            cds_seq += seq

        from Bio.Seq import Seq
        if cds_features[0].strand == '-':
            cds_seq = str(Seq(cds_seq).reverse_complement())

        if len(cds_seq) % 3 != 0:
            results['frameshift'] += 1
            continue

        protein = str(Seq(cds_seq).translate())
        if not protein.startswith('M'):
            results['no_start'] += 1
        elif not protein.endswith('*'):
            results['no_stop'] += 1
        elif protein[:-1].count('*') > 0:
            results['internal_stop'] += 1
        else:
            results['valid'] += 1

    print(f'\n=== ORF Validation ===')
    print(f'Total transcripts checked: {total}')
    for key, val in results.items():
        pct = val / total * 100 if total > 0 else 0
        print(f'  {key}: {val} ({pct:.1f}%)')

    return results


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 3:
        print('Usage: compare_annotations.py <reference.gff3> <transferred.gff3> [target.fasta]')
        sys.exit(1)

    ref_stats, tgt_stats = compare_annotations(sys.argv[1], sys.argv[2])

    if len(sys.argv) > 3:
        validate_orfs(sys.argv[2], sys.argv[3])
