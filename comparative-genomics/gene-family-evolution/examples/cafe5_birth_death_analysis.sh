#!/usr/bin/env bash
# Reference: CAFE5 5.1.0+, treePL 1.0+, Python 3.10+, pandas 2.2+, statsmodels 0.14+, R 4.4+ | Verify API if version differs
# CAFE5 birth-death analysis with annotation-error mode + FDR correction.

set -euo pipefail

ORTHOFINDER_OUT=${1:?usage: $0 ORTHOFINDER_OUT SPECIES_TREE_NWK OUTPUT_DIR}
SPECIES_TREE_NWK=${2:?missing species tree (Newick with substitution branches)}
OUTPUT_DIR=${3:?missing output dir}
THREADS=${THREADS:-16}

mkdir -p "$OUTPUT_DIR"/{cafe_input,tree,results}

echo "[1/5] Build CAFE5 input from OrthoFinder HOGs"
# v3 layout: Phylogenetic_Hierarchical_Orthogroups/N0.tsv
HOG_FILE=$(find "$ORTHOFINDER_OUT" -name 'N0.tsv' | head -1)
if [ -z "$HOG_FILE" ]; then
    echo "ERROR: HOG N0.tsv not found under $ORTHOFINDER_OUT"
    exit 1
fi

python3 - <<PY
import pandas as pd

hog = pd.read_csv('$HOG_FILE', sep='\t')
sp_cols = [c for c in hog.columns if c not in ('HOG', 'OG', 'Gene Tree Parent Clade')]

# Build count matrix: families x species
count_matrix = pd.DataFrame(index=hog['HOG'])
for sp in sp_cols:
    count_matrix[sp] = hog[sp].fillna('').apply(
        lambda x: len(str(x).split(',')) if x and ',' in str(x) else (1 if str(x).strip() else 0)
    )

# Filter to multi-copy families (max >= 2)
multi_copy = count_matrix[count_matrix.max(axis=1) >= 2]
print(f'Multi-copy families: {len(multi_copy)} / {len(count_matrix)}')

# CAFE5 format: Family ID + per-species counts
multi_copy.reset_index(inplace=True)
multi_copy.rename(columns={'index': 'family_id'}, inplace=True)
multi_copy.insert(0, 'Description', '(null)')

multi_copy.to_csv('$OUTPUT_DIR/cafe_input/family_counts.tsv', sep='\t', index=False)
print('CAFE5 input: $OUTPUT_DIR/cafe_input/family_counts.tsv')
PY

echo "[2/5] Time-calibrate species tree (ultrametric required for CAFE5)"
# Option 1: ape::chronos (R; default if no calibration)
Rscript - <<RS
library(ape)
tree <- read.tree('$SPECIES_TREE_NWK')
tree_u <- chronos(tree, lambda = 1, model = 'correlated')
write.tree(tree_u, '$OUTPUT_DIR/tree/tree_ultrametric.nwk')
RS

# Option 2 (alternative): treePL for fast NP-style calibration
# treepl config.txt

echo "[3/5] Run CAFE5 with gamma rate categories"
cafe5 \
    -i "$OUTPUT_DIR/cafe_input/family_counts.tsv" \
    -t "$OUTPUT_DIR/tree/tree_ultrametric.nwk" \
    -p \
    -k 4 \
    -c "$THREADS" \
    -o "$OUTPUT_DIR/results"

echo "[4/5] FDR correction across families"
python3 - <<PY
import pandas as pd
from statsmodels.stats.multitest import multipletests

# CAFE5 produces Base_family_results.txt with per-family p-values
results = pd.read_csv('$OUTPUT_DIR/results/Base_family_results.txt', sep='\t')
if 'pvalue' in results.columns:
    results['fdr'] = multipletests(results['pvalue'].fillna(1.0), method='fdr_bh')[1]
elif 'p-value' in results.columns:
    results['fdr'] = multipletests(results['p-value'].fillna(1.0), method='fdr_bh')[1]

# Identify significant expansions / contractions
significant = results[results['fdr'] < 0.05]
print(f'\nSignificant families (FDR < 0.05): {len(significant)} / {len(results)}')
significant.to_csv('$OUTPUT_DIR/results/significant_families_FDR.tsv', sep='\t', index=False)

# Top 20 by FDR
print(significant.sort_values('fdr').head(20))
PY

echo "[5/5] Per-clade lambda comparison (if model selection done)"
if [ -f "$OUTPUT_DIR/results/Base_clade_results.txt" ]; then
    cat "$OUTPUT_DIR/results/Base_clade_results.txt"
fi

echo "Done. Results: $OUTPUT_DIR/results"
