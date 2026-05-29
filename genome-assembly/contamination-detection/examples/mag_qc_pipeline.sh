#!/bin/bash
# Reference: CheckM2 1.0+, GUNC 1.0+, GTDB-Tk 2.4+, pandas 2.2+ | Verify API if version differs
# MAG/bin quality workflow: CheckM2 + GUNC TOGETHER + GTDB-Tk, gated on MIMAG.
# This is PROBLEM 2 (MAG quality). For PROBLEM 1 (foreign sequence in a single-organism
# assembly) run FCS-adaptor -> FCS-GX -> BlobToolKit instead - see fcs_screen.sh.
set -euo pipefail

BINS_DIR=$1
OUTPUT_DIR=$2
THREADS=${3:-16}
EXT=${4:-fa}

mkdir -p "$OUTPUT_DIR"

echo "=== CheckM2 (completeness + contamination via marker redundancy) ==="
checkm2 predict -i "$BINS_DIR" -o "$OUTPUT_DIR/checkm2" --threads "$THREADS" -x "$EXT"

echo "=== GUNC (chimerism via clade separation - catches disjoint-marker chimeras CheckM2 misses) ==="
gunc run -d "$BINS_DIR" -o "$OUTPUT_DIR/gunc" -t "$THREADS" -e ".$EXT"

echo "=== GTDB-Tk (taxonomy; DB release must match the binary) ==="
gtdbtk classify_wf --genome_dir "$BINS_DIR" --out_dir "$OUTPUT_DIR/gtdbtk" --extension "$EXT" --cpus "$THREADS"

echo "=== Merge CheckM2 x GUNC and gate on MIMAG ==="
OUTPUT_DIR="$OUTPUT_DIR" python3 << 'EOF'
import os
import pandas as pd

CONTAM_HQ = 5       # MIMAG high-quality: contamination < 5% (Bowers 2017); above this gene-content inference is unreliable
COMPLETE_HQ = 90    # MIMAG high-quality: completeness > 90%
CONTAM_MQ = 10      # MIMAG medium-quality: contamination < 10%
COMPLETE_MQ = 50    # MIMAG medium-quality: completeness >= 50%
RRS_TRUST = 0.5     # GUNC pass is trustworthy only when reference_representation_score > 0.5; below it 'pass' means 'can't tell'

out = os.environ['OUTPUT_DIR']
checkm = pd.read_csv(f'{out}/checkm2/quality_report.tsv', sep='\t')
gunc = pd.read_csv(f'{out}/gunc/GUNC.progenomes_2.1.maxCSS_level.tsv', sep='\t')
merged = checkm.merge(gunc, left_on='Name', right_on='genome', how='left')

merged['gunc_trustworthy'] = merged['reference_representation_score'] > RRS_TRUST
merged['mimag_high'] = ((merged['Completeness'] > COMPLETE_HQ) & (merged['Contamination'] < CONTAM_HQ) &
                        (merged['pass.GUNC'] == True) & merged['gunc_trustworthy'])
merged['mimag_medium'] = ((merged['Completeness'] >= COMPLETE_MQ) & (merged['Contamination'] < CONTAM_MQ) & ~merged['mimag_high'])
merged['chimera_suspect'] = (merged['Contamination'] < CONTAM_HQ) & (merged['pass.GUNC'] == False)

merged.to_csv(f'{out}/combined_qc.tsv', sep='\t', index=False)

print(f'MIMAG high-quality:   {int(merged["mimag_high"].sum())}')
print(f'MIMAG medium-quality: {int(merged["mimag_medium"].sum())}')
print(f'Chimera suspects (CheckM2-clean but GUNC-fail): {int(merged["chimera_suspect"].sum())}')
print(f'GUNC pass at low RRS (verdict unreliable): {int((~merged["gunc_trustworthy"] & (merged["pass.GUNC"] == True)).sum())}')
EOF

echo "=== Summary ==="
echo "Combined CheckM2 x GUNC table: $OUTPUT_DIR/combined_qc.tsv"
echo "Report the completeness/contamination PAIR with its GUNC verdict, never the contamination % alone."
