'''m6Anet workflow for ONT direct RNA m6A detection'''
import subprocess
from pathlib import Path
import pandas as pd

def run_m6anet_pipeline(fast5_dir, transcriptome_fa, output_dir, n_processes=8):
    '''Run complete m6Anet pipeline'''
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    data_dir = output_dir / 'dataprep'
    results_dir = output_dir / 'inference'

    # Step 1: Preprocess FAST5 files
    # Extracts signal features at DRACH motifs
    print('Running m6Anet dataprep...')
    subprocess.run([
        'm6anet', 'dataprep',
        '--input_dir', str(fast5_dir),
        '--output_dir', str(data_dir),
        '--reference', str(transcriptome_fa),
        '--n_processes', str(n_processes)
    ], check=True)

    # Step 2: Run m6A inference
    # Uses neural network to predict m6A probability
    print('Running m6Anet inference...')
    subprocess.run([
        'm6anet', 'inference',
        '--input_dir', str(data_dir),
        '--output_dir', str(results_dir),
        '--n_processes', str(n_processes)
    ], check=True)

    return results_dir

def filter_m6a_sites(results_dir, prob_threshold=0.9, min_coverage=20):
    '''Filter m6Anet results for high-confidence sites'''
    results_file = Path(results_dir) / 'data.site_proba.csv'
    df = pd.read_csv(results_file)

    # prob_threshold > 0.9: High confidence, few false positives
    # Lower (0.7-0.8) for more sensitive detection
    # min_coverage > 20: Require sufficient reads for reliable probability
    filtered = df[
        (df['probability_modified'] > prob_threshold) &
        (df['n_reads'] >= min_coverage)
    ]

    print(f'Total sites tested: {len(df)}')
    print(f'High-confidence m6A sites: {len(filtered)}')

    return filtered

def summarize_by_transcript(filtered_sites):
    '''Summarize m6A sites per transcript'''
    summary = filtered_sites.groupby('transcript_id').agg(
        n_m6a_sites=('position', 'count'),
        mean_probability=('probability_modified', 'mean'),
        total_coverage=('n_reads', 'sum')
    ).reset_index()

    return summary.sort_values('n_m6a_sites', ascending=False)

if __name__ == '__main__':
    import sys

    if len(sys.argv) < 4:
        print('Usage: m6anet_workflow.py <fast5_dir> <transcriptome.fa> <output_dir>')
        sys.exit(1)

    fast5_dir = sys.argv[1]
    transcriptome = sys.argv[2]
    output = sys.argv[3]

    results = run_m6anet_pipeline(fast5_dir, transcriptome, output)
    sites = filter_m6a_sites(results)
    sites.to_csv(Path(output) / 'm6a_high_confidence.csv', index=False)

    summary = summarize_by_transcript(sites)
    summary.to_csv(Path(output) / 'm6a_per_transcript.csv', index=False)
