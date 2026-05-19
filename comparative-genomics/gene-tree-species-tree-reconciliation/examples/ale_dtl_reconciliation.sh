#!/usr/bin/env bash
# Reference: ALE 1.0+, IQ-TREE 2.3.6+, Python 3.10+, ete4 4.1+, pandas 2.2+ | Verify API if version differs
# ALE probabilistic DTL reconciliation pipeline on orthogroup gene trees against a rooted species tree.
# Demonstrates: per-family UFBoot bootstrap, ALE encoding, ALEml_undated, per-branch event aggregation.

set -euo pipefail

ORTHOGROUP_DIR=${1:?usage: $0 ORTHOGROUP_DIR SPECIES_TREE OUTPUT_DIR}
SPECIES_TREE=${2:?missing species tree}
OUTPUT_DIR=${3:?missing output dir}
THREADS=${THREADS:-16}

mkdir -p "$OUTPUT_DIR"/{gene_trees,ale_encoded,reconciled}

echo "[1/4] Build UFBoot gene trees per orthogroup"
for og in "$ORTHOGROUP_DIR"/*.fa; do
    base=$(basename "$og" .fa)
    if [ ! -f "$OUTPUT_DIR/gene_trees/${base}.ufboot" ]; then
        iqtree2 -s "$og" -m TEST -B 1000 -nt 2 \
            --prefix "$OUTPUT_DIR/gene_trees/${base}"
    fi
done

echo "[2/4] Encode gene-tree distributions for ALE"
for ufb in "$OUTPUT_DIR/gene_trees"/*.ufboot; do
    cp "$ufb" "$OUTPUT_DIR/ale_encoded/"
    ALEobserve "$OUTPUT_DIR/ale_encoded/$(basename "$ufb")"
done

echo "[3/4] Reconcile against species tree"
for ale_file in "$OUTPUT_DIR/ale_encoded"/*.ale; do
    base=$(basename "$ale_file" .ale)
    out_uts="$OUTPUT_DIR/reconciled/${base}.uTs"
    if [ ! -f "$out_uts" ]; then
        ALEml_undated "$SPECIES_TREE" "$ale_file" \
            separators="|" \
            sample=100 \
            output_format=newick
        mv "${ale_file%.*}".*.uTs "$out_uts" 2>/dev/null || true
        mv "${ale_file%.*}".*.uml "$OUTPUT_DIR/reconciled/${base}.uml" 2>/dev/null || true
    fi
done

echo "[4/4] Aggregate per-branch events"
python3 - <<'PY'
import glob, os
from collections import defaultdict
import pandas as pd

reconciled_dir = os.environ.get('OUTPUT_DIR', '.') + '/reconciled'
agg = defaultdict(lambda: defaultdict(float))
for uts_path in glob.glob(f'{reconciled_dir}/*.uTs'):
    with open(uts_path) as fh:
        for ln in fh:
            if ln.startswith('#') or not ln.strip():
                continue
            parts = ln.split()
            if len(parts) >= 5:
                branch_id = parts[0]
                # branch  duplications  transfers  losses  originations
                agg[branch_id]['duplications'] += float(parts[1])
                agg[branch_id]['transfers']    += float(parts[2])
                agg[branch_id]['losses']       += float(parts[3])
                agg[branch_id]['originations'] += float(parts[4])

df = pd.DataFrame(agg).T.fillna(0).astype(int)
df.index.name = 'branch_id'
df.to_csv(f'{reconciled_dir}/../branch_events_summary.tsv', sep='\t')
print(df.sort_values('transfers', ascending=False).head(20))
PY

echo "Done. See $OUTPUT_DIR/branch_events_summary.tsv"
