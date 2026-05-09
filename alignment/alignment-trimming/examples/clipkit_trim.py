'''Trim alignment with ClipKIT and compute retention metrics.

ClipKIT modes (Steenwyk et al 2020 PLOS Bio):
- smart-gap (default): dynamic gap-fraction inflection point
- kpic-smart-gap: keep parsimony-informative + constant + smart-gap (recommended for trees)
- gappy: fixed gap threshold (-g 0.9 default)
- kpi: keep only parsimony-informative sites (most aggressive)
'''
# Reference: clipkit 2.1+ | Verify CLI flags if version differs

import subprocess
from Bio import AlignIO

def run_clipkit(input_fasta, output_fasta, mode='kpic-smart-gap', gap_threshold=None, log=True):
    cmd = ['clipkit', input_fasta, '-m', mode, '-o', output_fasta]
    if gap_threshold is not None:
        cmd += ['-g', str(gap_threshold)]
    if log:
        cmd += ['--log']
    subprocess.run(cmd, check=True)

def trimming_summary(input_fasta, output_fasta):
    original = AlignIO.read(input_fasta, 'fasta')
    trimmed = AlignIO.read(output_fasta, 'fasta')
    retention = trimmed.get_alignment_length() / original.get_alignment_length()
    return {
        'original_columns': original.get_alignment_length(),
        'trimmed_columns': trimmed.get_alignment_length(),
        'retention': retention,
        'removed': original.get_alignment_length() - trimmed.get_alignment_length(),
    }

if __name__ == '__main__':
    run_clipkit('input.fasta', 'trimmed.fasta', mode='kpic-smart-gap')
    summary = trimming_summary('input.fasta', 'trimmed.fasta')

    print(f'Original columns:  {summary["original_columns"]}')
    print(f'Trimmed columns:   {summary["trimmed_columns"]}')
    print(f'Retention:         {summary["retention"]*100:.1f}%')

    if summary['retention'] < 0.7:
        print('\nWARNING: Retention below 70%; trimming may degrade ML tree accuracy.')
        print('Tan, Muffato et al 2015 + Steenwyk 2020: aggressive trimming hurts on empirical data.')
        print('Consider switching to a less aggressive mode (smart-gap, kpic-gappy).')
