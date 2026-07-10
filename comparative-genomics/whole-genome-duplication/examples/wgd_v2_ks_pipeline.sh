#!/usr/bin/env bash
# Reference: wgd 2.0.31+, KsRates 1.1.3+, DupGen_finder 2024+, MCScanX 1.0+, PAML 4.10+, Python 3.10+ | Verify API if version differs
# wgd v2 + KsRates + DupGen_finder pipeline for WGD detection, dating, and duplication classification.

set -euo pipefail

# wgd v2 `dmd` takes CDS and translates internally (protein is optional and unused in a single-species
# paranome run), so this pipeline needs only the CDS FASTA -- do not require a protein FASTA it ignores.
CDS_FASTA=${1:?usage: $0 CDS_FASTA GFF OUTPUT_DIR}
GFF=${2:?missing GFF}
OUTPUT_DIR=${3:?missing output dir}
THREADS=${THREADS:-16}

mkdir -p "$OUTPUT_DIR"/{paranome,ks_dist,syn,mix,dupgen}

echo "[1/5] Build paranome (all-vs-all paralog identification)"
wgd dmd "$CDS_FASTA" -o "$OUTPUT_DIR/paranome" -t "$THREADS"

echo "[2/5] Compute Ks distribution"
wgd ksd "$OUTPUT_DIR/paranome/paranome.tsv" "$CDS_FASTA" \
    -o "$OUTPUT_DIR/ks_dist" \
    --aligner mafft --n-threads "$THREADS"
# Note: wgd ksd flags evolve across releases; verify with `wgd ksd --help`.

echo "[3/5] Synteny-anchored Ks (intra-genome)"
# Generate 4-col BED from GFF (chr  gene  start  end). Use Python for portability
# (the 3-arg match() form is gawk-only; BSD awk on macOS lacks it).
python3 -c "
import sys, re
with open('$GFF') as fh, open('$OUTPUT_DIR/syn/genome.bed', 'w') as out:
    for line in fh:
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) >= 9 and f[2] == 'gene':
            m = re.search(r'ID=([^;]+)', f[8])
            if m:
                out.write('\t'.join([f[0], m.group(1), f[3], f[4]]) + '\n')
"

wgd syn "$OUTPUT_DIR/paranome/paranome.tsv" "$OUTPUT_DIR/syn/genome.bed" "$CDS_FASTA" \
    -o "$OUTPUT_DIR/syn" \
    --min-block-size 5

echo "[4/5] Mixture model on Ks distribution"
python3 - <<PY
import pandas as pd
import numpy as np
from sklearn.mixture import GaussianMixture
import os

ksd_path = os.path.join('$OUTPUT_DIR', 'ks_dist', 'ks_distributions.tsv')
if not os.path.exists(ksd_path):
    ksd_path = next(iter([os.path.join('$OUTPUT_DIR', 'ks_dist', f) for f in os.listdir('$OUTPUT_DIR/ks_dist') if 'ks' in f.lower()]))

df = pd.read_csv(ksd_path, sep='\t')
# Filter saturated
ks_filt = df[(df['Ks'] > 0) & (df['Ks'] < 2)]['Ks'].values

# BIC-based model selection 1..5 components
bics = {}
for n in range(1, 6):
    gmm = GaussianMixture(n_components=n, random_state=42)
    gmm.fit(ks_filt.reshape(-1, 1))
    bics[n] = gmm.bic(ks_filt.reshape(-1, 1))

best_n = min(bics, key=bics.get)
print(f'BIC-selected components: {best_n}')

# Fit best model and report peak means
gmm = GaussianMixture(n_components=best_n, random_state=42).fit(ks_filt.reshape(-1, 1))
means = gmm.means_.flatten()
sorted_peaks = sorted(zip(means, gmm.weights_), key=lambda x: x[0])
print('Ks peak means (weights):')
for mean, weight in sorted_peaks:
    print(f'  Ks = {mean:.3f}, weight = {weight:.3f}')
PY

echo "[5/5] Classify duplications with DupGen_finder"
# Requires MCScanX collinearity already; if not, run:
# python -m jcvi.compara.catalog ortholog ... -> .collinearity file
# Then:
# perl DupGen_finder.pl -i <input> -t <species> -c <collinearity> -o $OUTPUT_DIR/dupgen
echo "Run DupGen_finder.pl manually for paralog classification:"
echo "  perl path/to/DupGen_finder.pl -i $OUTPUT_DIR/paranome -t species -c collinearity.col -o $OUTPUT_DIR/dupgen"

echo "Done. Results in $OUTPUT_DIR"
