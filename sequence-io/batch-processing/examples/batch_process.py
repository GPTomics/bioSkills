#!/usr/bin/env python3
'''Batch process multiple sequence files'''

from Bio import SeqIO
from pathlib import Path
import multiprocessing as mp
from functools import partial

def process_file(filepath, output_dir, min_length=100):
    '''Process a single FASTA/FASTQ file'''
    output_path = Path(output_dir) / filepath.name
    format_type = 'fastq' if filepath.suffix in ['.fq', '.fastq'] else 'fasta'

    count = 0
    with open(output_path, 'w') as out:
        for record in SeqIO.parse(filepath, format_type):
            if len(record.seq) >= min_length:
                SeqIO.write(record, out, format_type)
                count += 1

    return filepath.name, count

def batch_process(input_dir, output_dir, pattern='*.fasta', min_length=100, threads=4):
    '''Process multiple files in parallel'''
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    files = list(input_dir.glob(pattern))
    print(f'Processing {len(files)} files...')

    process_func = partial(process_file, output_dir=output_dir, min_length=min_length)

    with mp.Pool(threads) as pool:
        results = pool.map(process_func, files)

    for filename, count in results:
        print(f'  {filename}: {count} sequences')

    print(f'Output in: {output_dir}')

if __name__ == '__main__':
    import sys
    input_dir = sys.argv[1] if len(sys.argv) > 1 else 'input'
    output_dir = sys.argv[2] if len(sys.argv) > 2 else 'output'
    batch_process(input_dir, output_dir)
