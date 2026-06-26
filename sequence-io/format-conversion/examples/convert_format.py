#!/usr/bin/env python3
'''Convert between sequence file formats, flagging silent data loss'''
# Reference: biopython 1.83+ | Verify API if version differs

import sys
import gzip
from pathlib import Path
from Bio import SeqIO

FORMAT_MAP = {
    '.fasta': 'fasta', '.fa': 'fasta', '.fna': 'fasta',
    '.fastq': 'fastq', '.fq': 'fastq',
    '.gb': 'genbank', '.gbk': 'genbank',
    '.embl': 'embl',
    '.phy': 'phylip', '.phylip': 'phylip',
    '.aln': 'clustal',
    '.sth': 'stockholm',
}

LOSSY_TO_FASTA = {'genbank', 'embl', 'stockholm'}

def detect_format(filepath):
    '''Infer SeqIO format string from extension, ignoring a .gz suffix'''
    suffixes = [s for s in Path(filepath).suffixes if s != '.gz']
    if suffixes:
        return FORMAT_MAP.get(suffixes[-1].lower(), 'fasta')
    return 'fasta'

def smart_open(filepath, mode='r'):
    '''Open a file as a text handle, transparently handling gzip'''
    if str(filepath).endswith('.gz'):
        return gzip.open(filepath, mode + 't')
    return open(filepath, mode)

def warn_if_lossy(in_format, out_format):
    '''Print a warning when the target format cannot hold the source information'''
    if out_format == 'fasta' and in_format in LOSSY_TO_FASTA:
        print(f'WARNING: {in_format} -> fasta drops all features, qualifiers, annotations, and dbxrefs')
    if out_format == 'fasta' and in_format in ('fastq', 'fastq-sanger', 'fastq-illumina', 'fastq-solexa'):
        print('WARNING: fastq -> fasta drops per-base quality scores')

def convert(input_path, output_path, in_format=None, out_format=None):
    '''Convert a sequence file; uses SeqIO.convert when no gzip handle is needed'''
    if in_format is None:
        in_format = detect_format(input_path)
    if out_format is None:
        out_format = detect_format(output_path)
    warn_if_lossy(in_format, out_format)

    gz = str(input_path).endswith('.gz') or str(output_path).endswith('.gz')
    if gz:
        with smart_open(input_path, 'r') as fin, smart_open(output_path, 'w') as fout:
            count = SeqIO.write(SeqIO.parse(fin, in_format), fout, out_format)
    else:
        count = SeqIO.convert(input_path, in_format, output_path, out_format)

    print(f'Converted {count} records: {in_format} -> {out_format}')
    return count

def reencode_fastq(input_path, output_path, in_variant, out_variant):
    '''Re-encode FASTQ quality offsets between named variants (never auto-detected)'''
    count = SeqIO.convert(input_path, in_variant, output_path, out_variant)
    print(f'Re-encoded {count} records: {in_variant} -> {out_variant}')
    return count

if __name__ == '__main__':
    if len(sys.argv) >= 3:
        convert(sys.argv[1], sys.argv[2])
    else:
        print('Usage: convert_format.py input.gb output.fasta')
