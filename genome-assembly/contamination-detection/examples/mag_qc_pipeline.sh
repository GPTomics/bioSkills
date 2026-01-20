#!/bin/bash
set -euo pipefail

GENOMES_DIR=$1
OUTPUT_DIR=$2
THREADS=${3:-16}

mkdir -p "$OUTPUT_DIR"/{checkm2,gunc,gtdbtk,filtered}

echo "=== Running CheckM2 ==="
checkm2 predict \
    --input "$GENOMES_DIR" \
    --output-directory "$OUTPUT_DIR/checkm2" \
    --threads "$THREADS" \
    --extension fa

echo "=== Running GUNC ==="
gunc run \
    -d "$GENOMES_DIR" \
    -o "$OUTPUT_DIR/gunc" \
    -t "$THREADS" \
    -e .fa

echo "=== Running GTDB-Tk ==="
gtdbtk classify_wf \
    --genome_dir "$GENOMES_DIR" \
    --out_dir "$OUTPUT_DIR/gtdbtk" \
    --extension fa \
    --cpus "$THREADS"

echo "=== Filtering high-quality MAGs ==="
python3 << 'EOF'
import pandas as pd
import os
import shutil

output_dir = os.environ.get('OUTPUT_DIR', 'qc_output')
genomes_dir = os.environ.get('GENOMES_DIR', 'genomes')

checkm = pd.read_csv(f'{output_dir}/checkm2/quality_report.tsv', sep='\t')
gunc = pd.read_csv(f'{output_dir}/gunc/GUNC.progenomes_2.1.maxCSS_level.tsv', sep='\t')

merged = checkm.merge(gunc, left_on='Name', right_on='genome', how='left')
merged['pass.GUNC'] = merged['pass.GUNC'].fillna(True)

hq = merged[(merged['Completeness'] > 90) & (merged['Contamination'] < 5) & (merged['pass.GUNC'] == True)]
mq = merged[(merged['Completeness'] >= 50) & (merged['Contamination'] < 10)]

hq.to_csv(f'{output_dir}/high_quality_mags.tsv', sep='\t', index=False)
mq.to_csv(f'{output_dir}/medium_quality_mags.tsv', sep='\t', index=False)

os.makedirs(f'{output_dir}/filtered/high_quality', exist_ok=True)
os.makedirs(f'{output_dir}/filtered/medium_quality', exist_ok=True)

for genome in hq['Name']:
    src = f'{genomes_dir}/{genome}.fa'
    if os.path.exists(src):
        shutil.copy(src, f'{output_dir}/filtered/high_quality/')

print(f'High-quality MAGs: {len(hq)}')
print(f'Medium-quality MAGs: {len(mq)}')
EOF

echo "=== Summary ==="
echo "Results in: $OUTPUT_DIR"
echo "  - CheckM2 results: $OUTPUT_DIR/checkm2/quality_report.tsv"
echo "  - GUNC results: $OUTPUT_DIR/gunc/"
echo "  - GTDB-Tk results: $OUTPUT_DIR/gtdbtk/"
echo "  - Filtered MAGs: $OUTPUT_DIR/filtered/"
echo "  - Quality summary: $OUTPUT_DIR/high_quality_mags.tsv"
