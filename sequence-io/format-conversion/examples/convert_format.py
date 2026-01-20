#!/usr/bin/env python3
'''Convert between sequence file formats'''

from Bio import SeqIO
import gzip
from pathlib import Path

FORMAT_MAP = {
    '.fasta': 'fasta', '.fa': 'fasta', '.fna': 'fasta',
    '.fastq': 'fastq', '.fq': 'fastq',
    '.gb': 'genbank', '.gbk': 'genbank',
    '.embl': 'embl',
    '.phy': 'phylip', '.phylip': 'phylip',
    '.aln': 'clustal',
    '.sth': 'stockholm',
}

def detect_format(filepath):
    '''Detect format from extension'''
    path = Path(filepath)
    suffixes = [s for s in path.suffixes if s != '.gz']
    if suffixes:
        return FORMAT_MAP.get(suffixes[-1].lower(), 'fasta')
    return 'fasta'

def smart_open(filepath, mode='r'):
    '''Open file, handling gzip'''
    if str(filepath).endswith('.gz'):
        return gzip.open(filepath, mode + 't')
    return open(filepath, mode)

def convert(input_path, output_path, in_format=None, out_format=None):
    '''Convert between formats'''
    if in_format is None:
        in_format = detect_format(input_path)
    if out_format is None:
        out_format = detect_format(output_path)

    with smart_open(input_path, 'r') as fin, smart_open(output_path, 'w') as fout:
        records = SeqIO.parse(fin, in_format)
        count = SeqIO.write(records, fout, out_format)

    print(f'Converted {count} records: {in_format} -> {out_format}')
    return count

if __name__ == '__main__':
    import sys
    if len(sys.argv) >= 3:
        convert(sys.argv[1], sys.argv[2])
    else:
        print('Usage: convert_format.py input.fasta output.gb')
