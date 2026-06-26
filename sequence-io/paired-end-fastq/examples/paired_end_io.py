#!/usr/bin/env python3
'''Synchronized paired-end FASTQ filtering, interleaving, and deinterleaving.

Keeps R1/R2 mates in lockstep: a pair survives only if both mates pass, and a
read whose mate was discarded is routed to a separate orphan file rather than
silently desyncing the two streams (the #1 paired-end correctness trap).'''
# Reference: biopython 1.83+ | Verify API if version differs

import gzip
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def opener(path, mode):
    '''Open a plain or gzipped FASTQ in text mode based on the .gz suffix.'''
    if str(path).endswith('.gz'):
        return gzip.open(path, mode)
    return open(path, mode)

def mate_key(record):
    '''Shared mate ID: record.id drops anything after the first space (CASAVA 1.8+); rsplit drops a /1 or /2 suffix (pre-1.8).'''
    return record.id.rsplit('/', 1)[0]

def parse_paired(r1_path, r2_path):
    with opener(r1_path, 'rt') as f1, opener(r2_path, 'rt') as f2:
        for r1, r2 in zip(SeqIO.parse(f1, 'fastq'), SeqIO.parse(f2, 'fastq')):
            if mate_key(r1) != mate_key(r2):
                raise ValueError(f'Pair mismatch: {r1.id} vs {r2.id}')
            yield r1, r2

def mean_quality(record):
    return sum(record.letter_annotations['phred_quality']) / len(record.seq)

def filter_pairs_synced(r1_in, r2_in, r1_out, r2_out, r1_orphan, r2_orphan, min_qual=25):
    '''Keep a pair only if both mates pass; route lone survivors to orphan files.

    min_qual=25 is a common Phred cutoff (~0.3% error); tune per experiment.'''
    counts = {'paired': 0, 'r1_orphan': 0, 'r2_orphan': 0}
    with opener(r1_out, 'wt') as p1, opener(r2_out, 'wt') as p2, \
         opener(r1_orphan, 'wt') as o1, opener(r2_orphan, 'wt') as o2:
        for r1, r2 in parse_paired(r1_in, r2_in):
            r1_ok = mean_quality(r1) >= min_qual
            r2_ok = mean_quality(r2) >= min_qual
            if r1_ok and r2_ok:
                SeqIO.write(r1, p1, 'fastq')
                SeqIO.write(r2, p2, 'fastq')
                counts['paired'] += 1
            elif r1_ok:
                SeqIO.write(r1, o1, 'fastq')
                counts['r1_orphan'] += 1
            elif r2_ok:
                SeqIO.write(r2, o2, 'fastq')
                counts['r2_orphan'] += 1
    return counts

def interleave_pairs(r1_in, r2_in, out_path):
    def records():
        for r1, r2 in parse_paired(r1_in, r2_in):
            yield r1
            yield r2
    with opener(out_path, 'wt') as out:
        return SeqIO.write(records(), out, 'fastq') // 2

def deinterleave(in_path, r1_out, r2_out):
    pairs = 0
    with opener(in_path, 'rt') as fin, opener(r1_out, 'wt') as f1, opener(r2_out, 'wt') as f2:
        for i, record in enumerate(SeqIO.parse(fin, 'fastq')):
            if i % 2 == 0:
                SeqIO.write(record, f1, 'fastq')
            else:
                SeqIO.write(record, f2, 'fastq')
                pairs += 1
    return pairs

def _make_record(name, seq, qual):
    rec = SeqRecord(Seq(seq), id=name, description=f'{name} 1:N:0:ATCACG')
    rec.letter_annotations['phred_quality'] = qual
    return rec

def _demo():
    work = Path(__file__).parent / '_demo'
    work.mkdir(exist_ok=True)
    r1, r2 = work / 'R1.fastq', work / 'R2.fastq'

    high = [38] * 20
    low = [10] * 20  # well below min_qual to force an orphan
    with open(r1, 'w') as f1, open(r2, 'w') as f2:
        SeqIO.write(_make_record('read1', 'ACGT' * 5, high), f1, 'fastq')
        SeqIO.write(_make_record('read1', 'TTGG' * 5, high), f2, 'fastq')
        SeqIO.write(_make_record('read2', 'GGCC' * 5, high), f1, 'fastq')      # R1 passes
        SeqIO.write(_make_record('read2', 'AATT' * 5, low), f2, 'fastq')       # R2 fails -> R1 orphan

    counts = filter_pairs_synced(r1, r2, work / 'p_R1.fastq', work / 'p_R2.fastq',
                                 work / 'orphan_R1.fastq', work / 'orphan_R2.fastq')
    print(f'Filter: {counts}')

    n = interleave_pairs(work / 'p_R1.fastq', work / 'p_R2.fastq', work / 'inter.fastq')
    print(f'Interleaved {n} pairs')
    back = deinterleave(work / 'inter.fastq', work / 'back_R1.fastq', work / 'back_R2.fastq')
    print(f'Deinterleaved {back} pairs')

    for child in work.iterdir():
        child.unlink()
    work.rmdir()

if __name__ == '__main__':
    _demo()
