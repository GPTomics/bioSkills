#!/usr/bin/env python3
# Reference: deepTools 3.5+, pysam 0.22+, pandas 2.2+, PreSeq 3.2+ | Verify with `deeptools --version`, `pysam.__version__`, `preseq --version` if installed releases differ.
# MeRIP-seq QC: replicate Spearman concordance, plotFingerprint JS distances, saturation curves via PreSeq.
# Inputs: sorted indexed BAM files for IP and Input libraries in aligned/.

import subprocess
from pathlib import Path

import pandas as pd

ALIGNED_DIR = Path('aligned')
QC_DIR = Path('qc')
COMPLEXITY_DIR = Path('complexity')
TRACK_DIR = Path('tracks')
BIN_SIZE_CORRELATION = 10000
BIN_SIZE_TRACKS = 25
PSEUDOCOUNT = 1
THREADS = 8

QC_DIR.mkdir(exist_ok=True)
COMPLEXITY_DIR.mkdir(exist_ok=True)
TRACK_DIR.mkdir(exist_ok=True)

ip_bams = sorted(ALIGNED_DIR.glob('IP_*_Aligned.sortedByCoord.out.bam'))
input_bams = sorted(ALIGNED_DIR.glob('Input_*_Aligned.sortedByCoord.out.bam'))
all_bams = ip_bams + input_bams
labels = [b.name.split('_Aligned')[0] for b in all_bams]

# Replicate concordance: 10 kb bins -> Spearman heatmap.
subprocess.run(['multiBamSummary', 'bins',
    '--bamfiles', *map(str, all_bams),
    '--binSize', str(BIN_SIZE_CORRELATION),
    '--numberOfProcessors', str(THREADS),
    '--outRawCounts', str(QC_DIR / 'raw_bin_counts.tab'),
    '-o', str(QC_DIR / 'cov.npz')], check=True)

subprocess.run(['plotCorrelation',
    '--corData', str(QC_DIR / 'cov.npz'),
    '--corMethod', 'spearman',
    '--skipZeros',
    '--whatToPlot', 'heatmap',
    '--colorMap', 'RdYlBu_r',
    '--plotNumbers',
    '-o', str(QC_DIR / 'replicate_correlation.pdf')], check=True)

# IP enrichment QC via plotFingerprint.
subprocess.run(['plotFingerprint',
    '--bamfiles', *map(str, all_bams),
    '--labels', *labels,
    '--numberOfProcessors', str(THREADS),
    '--skipZeros',
    '--outQualityMetrics', str(QC_DIR / 'fingerprint_metrics.tab'),
    '-o', str(QC_DIR / 'fingerprint.pdf')], check=True)

# Library complexity / saturation per library. Subsampling for cross-condition comparison drives off these curves.
for bam in all_bams:
    name = bam.stem.split('_Aligned')[0]
    subprocess.run(['preseq', 'c_curve', '-B', '-o', str(COMPLEXITY_DIR / f'{name}_c_curve.txt'), str(bam)], check=True)
    subprocess.run(['preseq', 'lc_extrap', '-B', '-o', str(COMPLEXITY_DIR / f'{name}_lc_extrap.txt'), str(bam)], check=True)

# IP-over-Input log2 bigWig tracks for each matched IP/Input pair.
pairs = list(zip(ip_bams, input_bams))
for ip, inp in pairs:
    out = TRACK_DIR / f'{ip.stem.split("_Aligned")[0]}_over_{inp.stem.split("_Aligned")[0]}.bw'
    subprocess.run(['bamCompare',
        '-b1', str(ip),
        '-b2', str(inp),
        '--operation', 'log2',
        '--pseudocount', str(PSEUDOCOUNT),
        '--binSize', str(BIN_SIZE_TRACKS),
        '--normalizeUsing', 'CPM',
        '--numberOfProcessors', str(THREADS),
        '-o', str(out)], check=True)

# Summarise fingerprint metrics; flag failed IPs (low IP-vs-input JS distance).
fp = pd.read_csv(QC_DIR / 'fingerprint_metrics.tab', sep='\t')
ip_rows = fp[fp['Sample'].str.startswith('IP_')].copy()
ip_rows['flag'] = ip_rows['JS Distance'].apply(lambda j: 'POSSIBLE_FAILED_IP' if j < 0.5 else 'OK')
ip_rows[['Sample', 'JS Distance', 'flag']].to_csv(QC_DIR / 'ip_qc_summary.tsv', sep='\t', index=False)

print('QC complete. Inspect:')
print(f'  - {QC_DIR / "replicate_correlation.pdf"} (Spearman matrix; IP-IP within condition >= 0.85)')
print(f'  - {QC_DIR / "fingerprint.pdf"} (Lorenz curves; IP should sit below diagonal)')
print(f'  - {COMPLEXITY_DIR}/*lc_extrap.txt (saturation curves; rarefy to common depth before peak calling)')
print(f'  - {QC_DIR / "ip_qc_summary.tsv"} (failed-IP flags)')
