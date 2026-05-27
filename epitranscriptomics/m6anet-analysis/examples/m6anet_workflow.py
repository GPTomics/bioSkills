#!/usr/bin/env python3
# Reference: m6anet 2.1+, nanopolish 0.14+, minimap2 2.26+, samtools 1.19+, pandas 2.2+ | Verify with `pip show m6anet`, `nanopolish --version`, `minimap2 --version`, `samtools --version`, `import pandas; pandas.__version__` if installed releases differ.
# Full m6Anet pipeline: basecalled FASTQ -> transcriptome alignment -> nanopolish eventalign -> m6anet dataprep -> m6anet inference.
# CRITICAL: minimap2 to TRANSCRIPTOME (not genome); nanopolish --scale-events --signal-index are the m6Anet-required flags (--samples --print-read-names are needed only by other downstream tools).

import subprocess
from pathlib import Path

import pandas as pd

READS_FASTQ = Path('reads.fastq')
POD5_DIR = Path('pod5')
TRANSCRIPTOME_FA = Path('refs/transcriptome.fa')
ALIGNED_BAM = Path('aligned.bam')
EVENTALIGN_OUT = Path('eventalign.txt')
DATAPREP_DIR = Path('m6anet_data')
INFERENCE_DIR = Path('m6anet_results')

THREADS = 12
PROBABILITY_THRESHOLD = 0.9
MIN_COVERAGE = 20

DATAPREP_DIR.mkdir(exist_ok=True)
INFERENCE_DIR.mkdir(exist_ok=True)

# Step 1: minimap2 transcriptome alignment with ONT DRS flags. -uf forces forward strand (DRS is directional); -k14 is recommended k-mer.
with open('aligned.sam', 'w') as out:
    subprocess.run(['minimap2',
        '-ax', 'map-ont',
        '-uf',
        '-k14',
        '--secondary=no',
        '-t', str(THREADS),
        str(TRANSCRIPTOME_FA),
        str(READS_FASTQ)], check=True, stdout=out)

subprocess.run(['samtools', 'sort', '-@', '8', '-o', str(ALIGNED_BAM), 'aligned.sam'], check=True)
subprocess.run(['samtools', 'index', str(ALIGNED_BAM)], check=True)
Path('aligned.sam').unlink()

# Step 2: nanopolish index + eventalign. Required flags pinned to m6Anet trained input format.
subprocess.run(['nanopolish', 'index', '-d', str(POD5_DIR), str(READS_FASTQ)], check=True)

with open(EVENTALIGN_OUT, 'w') as out:
    subprocess.run(['nanopolish', 'eventalign',
        '--reads', str(READS_FASTQ),
        '--bam', str(ALIGNED_BAM),
        '--genome', str(TRANSCRIPTOME_FA),
        '--scale-events',
        '--signal-index',
        '--threads', str(THREADS),
        '--summary', 'nanopolish_summary.tsv'], check=True, stdout=out)

# Step 3: m6anet dataprep + inference (v2+ subcommand syntax).
subprocess.run(['m6anet', 'dataprep',
    '--eventalign', str(EVENTALIGN_OUT),
    '--out_dir', str(DATAPREP_DIR),
    '--n_processes', str(THREADS)], check=True)

subprocess.run(['m6anet', 'inference',
    '--input_dir', str(DATAPREP_DIR),
    '--out_dir', str(INFERENCE_DIR),
    '--n_processes', '4',
    '--num_iterations', '1000'], check=True)

# Step 4: Filter sites and aggregate per-transcript.
sites = pd.read_csv(INFERENCE_DIR / 'data.site_proba.csv')

print(f'Total DRACH sites tested: {len(sites)}')

filtered = sites[(sites['n_reads'] >= MIN_COVERAGE) &
                 (sites['probability_modified'] >= PROBABILITY_THRESHOLD)]

print(f'High-confidence m6A sites (n_reads >= {MIN_COVERAGE}, prob >= {PROBABILITY_THRESHOLD}): {len(filtered)}')

per_transcript = (filtered
    .groupby('transcript_id')
    .agg(n_high_conf_sites=('transcript_position', 'count'),
         mean_probability=('probability_modified', 'mean'),
         mean_mod_ratio=('mod_ratio', 'mean'),
         total_coverage=('n_reads', 'sum'))
    .reset_index()
    .sort_values('n_high_conf_sites', ascending=False))

filtered.to_csv(INFERENCE_DIR / 'high_confidence_sites.tsv', sep='\t', index=False)
per_transcript.to_csv(INFERENCE_DIR / 'm6a_per_transcript.tsv', sep='\t', index=False)

print('Outputs:')
print(f'  - {INFERENCE_DIR / "data.site_proba.csv"} (all per-site predictions)')
print(f'  - {INFERENCE_DIR / "high_confidence_sites.tsv"} (filtered)')
print(f'  - {INFERENCE_DIR / "m6a_per_transcript.tsv"} (top modified transcripts)')
print('Next: cross-validate top hits against GLORI / m6A-SAC-seq for absolute stoichiometry.')
