'''Trim alignment with BMGE (entropy + BLOSUM62 context).

BMGE (Criscuolo & Gribaldo 2010 BMC Evol Biol) is the standard trimmer in deep
prokaryotic phylogenomics (GToTree pipeline). Recommended `-h 0.4 -g 0.2` for
deep phylogenies; `-h 0.6` retains more sites for shallower datasets.
'''
# Reference: BMGE 1.12+ | Verify CLI flags if version differs

import subprocess
from Bio import AlignIO

def run_bmge(input_fasta, output_fasta, sequence_type='AA', entropy=0.5, gap_threshold=0.2, matrix=None):
    cmd = [
        'java', '-Xmx8g', '-jar', 'BMGE.jar',
        '-i', input_fasta,
        '-of', output_fasta,
        '-t', sequence_type,
        '-h', str(entropy),
        '-g', str(gap_threshold),
    ]
    if matrix:
        cmd += ['-m', matrix]
    subprocess.run(cmd, check=True)

def trimming_summary(input_fasta, output_fasta):
    original = AlignIO.read(input_fasta, 'fasta')
    trimmed = AlignIO.read(output_fasta, 'fasta')
    return {
        'original': original.get_alignment_length(),
        'trimmed': trimmed.get_alignment_length(),
        'retention': trimmed.get_alignment_length() / original.get_alignment_length(),
    }

if __name__ == '__main__':
    run_bmge('input.fasta', 'trimmed.fasta', sequence_type='AA', entropy=0.4, gap_threshold=0.2)
    summary = trimming_summary('input.fasta', 'trimmed.fasta')
    print('BMGE -h 0.4 -g 0.2 (deep prokaryotic phylogenomics):')
    print(f'  {summary["original"]} -> {summary["trimmed"]} columns ({summary["retention"]*100:.1f}% retained)')

    run_bmge('input.fasta', 'trimmed_shallow.fasta', sequence_type='AA', entropy=0.6, gap_threshold=0.2)
    summary_shallow = trimming_summary('input.fasta', 'trimmed_shallow.fasta')
    print('\nBMGE -h 0.6 (shallower phylogeny, more sites retained):')
    print(f'  {summary_shallow["original"]} -> {summary_shallow["trimmed"]} columns ({summary_shallow["retention"]*100:.1f}% retained)')
