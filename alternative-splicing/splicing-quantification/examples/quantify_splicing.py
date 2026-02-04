#!/usr/bin/env python3
'''
Splicing quantification using SUPPA2 and rMATS-turbo.
Demonstrates PSI calculation from RNA-seq data.
'''

import subprocess
import pandas as pd
from pathlib import Path


def run_suppa2_quantification(gtf_file, tpm_file, output_prefix, event_types=None):
    '''
    Quantify splicing using SUPPA2 from transcript TPM.

    Args:
        gtf_file: Path to GTF annotation
        tpm_file: Tab-separated TPM file (transcripts x samples)
        output_prefix: Prefix for output files
        event_types: List of event types (default: all)
    '''
    if event_types is None:
        event_types = ['SE', 'SS', 'MX', 'RI', 'FL']

    # Step 1: Generate splicing events from annotation
    subprocess.run([
        'suppa.py', 'generateEvents',
        '-i', gtf_file,
        '-o', output_prefix,
        '-f', 'ioe',
        '-e'] + event_types,
        check=True
    )

    # Step 2: Calculate PSI for each event type
    event_codes = {'SE': 'SE', 'SS': 'A5', 'MX': 'MX', 'RI': 'RI', 'FL': 'FL'}
    psi_files = {}

    for et in event_types:
        code = event_codes.get(et, et)
        ioe_file = f'{output_prefix}_{code}_strict.ioe'
        psi_output = f'{output_prefix}_psi_{code}'

        if Path(ioe_file).exists():
            subprocess.run([
                'suppa.py', 'psiPerEvent',
                '-i', ioe_file,
                '-e', tpm_file,
                '-o', psi_output
            ], check=True)
            psi_files[code] = f'{psi_output}.psi'

    return psi_files


def filter_reliable_events(psi_file, min_samples_with_coverage=0.5):
    '''
    Filter PSI matrix for reliable events.

    Args:
        psi_file: Path to PSI file from SUPPA2
        min_samples_with_coverage: Minimum fraction of samples with non-NA PSI
    '''
    psi = pd.read_csv(psi_file, sep='\t', index_col=0)

    # Remove events with too many missing values
    # NA indicates insufficient reads for PSI calculation
    coverage_frac = psi.notna().mean(axis=1)
    reliable = psi[coverage_frac >= min_samples_with_coverage]

    # Remove constitutive events (PSI always near 0 or 1)
    # PSI 0.1-0.9 range indicates true alternative splicing
    mean_psi = reliable.mean(axis=1)
    variable = reliable[(mean_psi > 0.1) & (mean_psi < 0.9)]

    print(f'Total events: {len(psi)}')
    print(f'After coverage filter: {len(reliable)}')
    print(f'After variability filter: {len(variable)}')

    return variable


def parse_rmats_output(rmats_dir, event_type='SE', min_junction_reads=20):
    '''
    Parse rMATS output and calculate mean PSI.

    Args:
        rmats_dir: Directory containing rMATS output
        event_type: SE, A5SS, A3SS, MXE, or RI
        min_junction_reads: Minimum total junction reads for reliable PSI
    '''
    jc_file = Path(rmats_dir) / f'{event_type}.MATS.JC.txt'
    df = pd.read_csv(jc_file, sep='\t')

    # IncLevel columns contain PSI values (comma-separated per replicate)
    inc_cols = [c for c in df.columns if c.startswith('IncLevel')]

    # Parse comma-separated values and calculate mean
    def parse_inc_levels(row):
        values = []
        for col in inc_cols:
            if pd.notna(row[col]):
                values.extend([float(x) for x in str(row[col]).split(',') if x != 'NA'])
        return values

    df['psi_values'] = df.apply(parse_inc_levels, axis=1)
    df['mean_PSI'] = df['psi_values'].apply(lambda x: sum(x) / len(x) if x else None)

    # Filter by junction read coverage
    # Sum inclusion and skipping junction counts
    df['total_reads'] = (
        df['IJC_SAMPLE_1'].str.split(',').apply(lambda x: sum(int(i) for i in x)) +
        df['SJC_SAMPLE_1'].str.split(',').apply(lambda x: sum(int(i) for i in x))
    )
    reliable = df[df['total_reads'] >= min_junction_reads].copy()

    print(f'Total {event_type} events: {len(df)}')
    print(f'Reliable events (>={min_junction_reads} reads): {len(reliable)}')

    return reliable[['ID', 'GeneID', 'geneSymbol', 'chr', 'strand',
                     'exonStart_0base', 'exonEnd', 'mean_PSI', 'total_reads']]


if __name__ == '__main__':
    # Example usage with SUPPA2
    # Requires: annotation.gtf, transcript_tpm.tsv from Salmon/kallisto
    gtf = 'annotation.gtf'
    tpm = 'transcript_tpm.tsv'

    # Run SUPPA2 quantification
    # psi_files = run_suppa2_quantification(gtf, tpm, 'splicing_events')

    # Filter for reliable events
    # reliable_se = filter_reliable_events('splicing_events_psi_SE.psi')
    # print(reliable_se.head())

    # Example with rMATS output
    # reliable_se = parse_rmats_output('rmats_output/', 'SE', min_junction_reads=20)
    # print(reliable_se.head())

    print('Run with actual data files to quantify splicing events')
