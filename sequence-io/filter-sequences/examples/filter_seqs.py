#!/usr/bin/env python3
'''Filter sequences by various criteria'''

from Bio import SeqIO
import gzip

def filter_by_length(records, min_length=0, max_length=float('inf')):
    '''Filter by sequence length'''
    for record in records:
        if min_length <= len(record) <= max_length:
            yield record

def filter_by_quality(records, min_mean_quality=20):
    '''Filter FASTQ by mean quality'''
    for record in records:
        quals = record.letter_annotations.get('phred_quality', [])
        if quals and sum(quals) / len(quals) >= min_mean_quality:
            yield record

def filter_by_gc(records, min_gc=0, max_gc=100):
    '''Filter by GC content'''
    for record in records:
        gc = (record.seq.count('G') + record.seq.count('C')) / len(record) * 100
        if min_gc <= gc <= max_gc:
            yield record

def filter_by_id(records, ids_to_keep):
    '''Filter by sequence ID'''
    ids_set = set(ids_to_keep)
    for record in records:
        if record.id in ids_set:
            yield record

def filter_fastq(input_path, output_path, min_length=50, min_quality=20):
    '''Filter FASTQ file'''
    opener_in = gzip.open if input_path.endswith('.gz') else open
    opener_out = gzip.open if output_path.endswith('.gz') else open

    with opener_in(input_path, 'rt') as fin, opener_out(output_path, 'wt') as fout:
        records = SeqIO.parse(fin, 'fastq')
        filtered = filter_by_length(records, min_length=min_length)
        filtered = filter_by_quality(filtered, min_mean_quality=min_quality)
        count = SeqIO.write(filtered, fout, 'fastq')

    print(f'Wrote {count} sequences to {output_path}')
    return count

if __name__ == '__main__':
    import sys
    input_path = sys.argv[1] if len(sys.argv) > 1 else 'input.fastq.gz'
    output_path = sys.argv[2] if len(sys.argv) > 2 else 'filtered.fastq.gz'
    filter_fastq(input_path, output_path)
