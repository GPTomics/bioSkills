'''Run MAFFT, MUSCLE5, or ClustalOmega based on dataset size, then summarize the alignment.'''
# Reference: MAFFT 7.520+, MUSCLE 5.1+, ClustalOmega 1.2.4+ | Verify CLI flags if version differs

import subprocess
from pathlib import Path
from Bio import AlignIO


def run_mafft(input_fasta, output_fasta, algorithm='linsi', threads=4):
    algo_flags = {
        'linsi': ['--localpair', '--maxiterate', '1000'],
        'ginsi': ['--globalpair', '--maxiterate', '1000'],
        'einsi': ['--genafpair', '--maxiterate', '1000'],
        'fftns2': ['--retree', '2'],
        'auto': ['--auto'],
    }
    cmd = ['mafft', '--thread', str(threads)] + algo_flags[algorithm] + [str(input_fasta)]
    with open(output_fasta, 'w') as out:
        subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, check=True)


def run_muscle(input_fasta, output_fasta, threads=4, super5=False):
    if super5:
        cmd = ['muscle', '-super5', str(input_fasta), '-output', str(output_fasta), '-threads', str(threads)]
    else:
        cmd = ['muscle', '-align', str(input_fasta), '-output', str(output_fasta), '-threads', str(threads)]
    subprocess.run(cmd, check=True)


def run_clustalo(input_fasta, output_fasta, threads=4, iterations=0):
    cmd = ['clustalo', '-i', str(input_fasta), '-o', str(output_fasta), '--threads=%d' % threads, '--force']
    if iterations > 0:
        cmd.extend(['--iter=%d' % iterations])
    subprocess.run(cmd, check=True)


def select_and_run(input_fasta, output_fasta, num_seqs=None, threads=4):
    if num_seqs is None:
        with open(input_fasta) as f:
            num_seqs = sum(1 for line in f if line.startswith('>'))

    if num_seqs <= 200:
        print(f'{num_seqs} sequences: using MAFFT L-INS-i (highest accuracy)')
        run_mafft(input_fasta, output_fasta, algorithm='linsi', threads=threads)
    elif num_seqs <= 10000:
        print(f'{num_seqs} sequences: using MAFFT FFT-NS-2 (good balance)')
        run_mafft(input_fasta, output_fasta, algorithm='fftns2', threads=threads)
    else:
        print(f'{num_seqs} sequences: using ClustalOmega (scales best)')
        run_clustalo(input_fasta, output_fasta, threads=threads)


def summarize_alignment(alignment_file, fmt='fasta'):
    aln = AlignIO.read(alignment_file, fmt)
    n_seqs = len(aln)
    n_cols = aln.get_alignment_length()
    total_gaps = sum(str(r.seq).count('-') for r in aln)
    total_positions = n_seqs * n_cols
    gap_fraction = total_gaps / total_positions

    gap_free_cols = sum(1 for i in range(n_cols) if '-' not in aln[:, i])

    print(f'Sequences: {n_seqs}')
    print(f'Alignment length: {n_cols} columns')
    print(f'Gap-free columns: {gap_free_cols} ({gap_free_cols/n_cols*100:.1f}%)')
    print(f'Overall gap fraction: {gap_fraction*100:.1f}%')


if __name__ == '__main__':
    input_file = 'sequences.fasta'
    output_file = 'aligned.fasta'

    select_and_run(input_file, output_file)
    summarize_alignment(output_file)
