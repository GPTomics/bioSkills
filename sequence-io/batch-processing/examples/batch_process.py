#!/usr/bin/env python3
'''Batch count, split, and index sequence files with memory-safe streaming'''
# Reference: biopython 1.83+ | Verify API if version differs

import tempfile
from pathlib import Path
from itertools import islice
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def count_streaming(filepath, format):
    '''Count records without materializing them (sum over a generator, not list()).'''
    return sum(1 for _ in SeqIO.parse(filepath, format))


def split_by_count(input_file, format, records_per_file, output_prefix):
    '''Split a file into chunks of records_per_file each, streaming via islice.'''
    records = SeqIO.parse(input_file, format)
    written = []
    file_num = 1
    while True:
        batch = list(islice(records, records_per_file))
        if not batch:
            break
        out = f'{output_prefix}_{file_num}.{format}'
        SeqIO.write(batch, out, format)
        written.append(out)
        file_num += 1
    return written


def make_demo_fastas(directory, n_files, seqs_per_file):
    '''Write small demo FASTA files so the example runs standalone.'''
    paths = []
    for fi in range(n_files):
        path = Path(directory) / f'sample{fi}.fasta'
        records = [SeqRecord(Seq('ACGT' * (10 + i)), id=f's{fi}_{i}', description='') for i in range(seqs_per_file)]
        SeqIO.write(records, path, 'fasta')
        paths.append(str(path))
    return paths


def main():
    with tempfile.TemporaryDirectory() as workdir:
        files = make_demo_fastas(workdir, n_files=3, seqs_per_file=5)

        for f in files:
            print(f'{Path(f).name}: {count_streaming(f, "fasta")} sequences')

        chunks = split_by_count(files[0], 'fasta', 2, str(Path(workdir) / 'chunk'))
        print(f'Split {Path(files[0]).name} into {len(chunks)} chunks of <=2 records')

        index_path = str(Path(workdir) / 'combined.idx')
        records = SeqIO.index_db(index_path, files, 'fasta')
        print(f'Indexed {len(records)} records across {len(files)} files')
        print(f'Random lookup s2_3 length: {len(records["s2_3"].seq)}')
        records.close()


if __name__ == '__main__':
    main()
