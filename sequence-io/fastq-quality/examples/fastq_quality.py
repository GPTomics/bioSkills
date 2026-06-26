#!/usr/bin/env python3
'''Read FASTQ quality scores and demonstrate the encoding/offset distinction'''
# Reference: biopython 1.83+ | Verify API if version differs

import gzip
import io
from collections import defaultdict
from Bio import SeqIO


def parse_fastq(filepath, fmt='fastq'):
    '''Stream FASTQ records, transparently handling gzip. fmt names the encoding explicitly.'''
    opener = gzip.open if str(filepath).endswith('.gz') else open
    with opener(filepath, 'rt') as handle:
        for record in SeqIO.parse(handle, fmt):
            yield record


def quality_stats(records):
    '''Per-read mean quality, length, and Q30 fraction from letter_annotations.'''
    means, lengths = [], []
    for record in records:
        quals = record.letter_annotations['phred_quality']
        means.append(sum(quals) / len(quals))
        lengths.append(len(record))
    n = len(means)
    return {'n_reads': n, 'mean_quality': sum(means) / n, 'mean_length': sum(lengths) / n,
            'q30_read_fraction': sum(m >= 30 for m in means) / n}


def per_position_mean(records):
    '''Mean quality at each read position across all reads.'''
    position_quals = defaultdict(list)
    for record in records:
        for pos, q in enumerate(record.letter_annotations['phred_quality']):
            position_quals[pos].append(q)
    return {pos: sum(qs) / len(qs) for pos, qs in sorted(position_quals.items())}


# A single FASTQ record whose quality line uses only ASCII >= 64 (the overlap region).
# The SAME bytes decode to different Phred scores under different offsets, which is why
# the encoding must come from instrument metadata and never be guessed.
overlap_fastq = '@read1\nACGTACGT\n+\nIIIIDDDD\n'


def encoding_distinction():
    '''Same bytes, two offsets: Phred+33 vs Phred+64 differ by exactly 31 (silent corruption).'''
    sanger = next(SeqIO.parse(io.StringIO(overlap_fastq), 'fastq'))
    illumina = next(SeqIO.parse(io.StringIO(overlap_fastq), 'fastq-illumina'))
    return sanger.letter_annotations['phred_quality'], illumina.letter_annotations['phred_quality']


# Solexa uses an ODDS score that can go negative (floor -5); it lives under 'solexa_quality'.
solexa_fastq = '@read1\nACGT\n+\n;;;;\n'


def solexa_scores():
    record = next(SeqIO.parse(io.StringIO(solexa_fastq), 'fastq-solexa'))
    return record.letter_annotations['solexa_quality']


if __name__ == '__main__':
    phred33, phred64 = encoding_distinction()
    print('Same quality bytes "IIIIDDDD" decoded two ways:')
    print(f'  as Phred+33 (fastq):          {phred33}')
    print(f'  as Phred+64 (fastq-illumina): {phred64}')
    print(f'  every score differs by exactly {phred33[0] - phred64[0]} -> silent 31-shift if mislabeled')
    print(f'\nSolexa ";;;;" decoded as fastq-solexa: {solexa_scores()} (odds score, negative allowed)')

    records = list(SeqIO.parse(io.StringIO(overlap_fastq * 3), 'fastq'))
    print(f'\nQuality stats: {quality_stats(records)}')
    print(f'Per-position mean: {per_position_mean(records)}')
