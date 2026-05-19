#!/usr/bin/env bash
# Reference: Dsuite 0.5+, TreeMix 1.13+, bcftools 1.21+, OptM 0.1.7+, R 4.4+, Python 3.10+, pandas 2.2+ | Verify API if version differs
# Dsuite ABBA-BABA + Fbranch + TreeMix introgression detection pipeline.

set -euo pipefail

VCF=${1:?usage: $0 VCF SETS_TSV SPECIES_TREE OUTPUT_DIR}
SETS_TSV=${2:?missing SETS.tsv (sample_id  population)}
SPECIES_TREE=${3:?missing species_tree.nwk}
OUTPUT_DIR=${4:?missing output dir}
THREADS=${THREADS:-16}

mkdir -p "$OUTPUT_DIR"/{dsuite,treemix}

echo "[1/4] Dsuite Dtrios - compute D-statistic for all population trios"
# Verify flags with `Dsuite Dtrios --help`:
#   -t/--tree species tree; -o/--out-prefix output prefix
#   -k/--no-f4-ratio skip f4-ratio; -c/--no-combine skip _combine.txt
# Fbranch is its own subcommand below (step 2).
Dsuite Dtrios \
    --tree="$SPECIES_TREE" \
    -o "$OUTPUT_DIR/dsuite/trios_run" \
    "$VCF" \
    "$SETS_TSV"

SETS_BASE=$(basename "$SETS_TSV" .tsv)

echo "[2/4] Dsuite Fbranch - tree-aware admixture mapping"
TREE_FILE="$OUTPUT_DIR/dsuite/trios_run_tree.txt"
if [ -f "$TREE_FILE" ]; then
    Dsuite Fbranch "$SPECIES_TREE" "$TREE_FILE" \
        > "$OUTPUT_DIR/dsuite/trios_run_fbranch.txt"
fi

echo "[3/4] Convert VCF to TreeMix input"
# Create population map
awk '{print $2, $1}' "$SETS_TSV" > "$OUTPUT_DIR/treemix/popmap.tsv"

# Use plink for allele-frequency extraction
plink --vcf "$VCF" \
    --make-bed \
    --out "$OUTPUT_DIR/treemix/plink" \
    --double-id --const-fid 0

plink --bfile "$OUTPUT_DIR/treemix/plink" \
    --freq --within "$OUTPUT_DIR/treemix/popmap.tsv" \
    --out "$OUTPUT_DIR/treemix/plink_freq" \
    --double-id

# Convert to TreeMix format
python3 - <<PY
import gzip
import pandas as pd

# Read plink .frq.strat
df = pd.read_csv('$OUTPUT_DIR/treemix/plink_freq.frq.strat', sep='\s+')
# Convert to TreeMix format
pivot = df.pivot_table(index='SNP', columns='CLST', values='MAC').fillna(0)
mac_total = df.pivot_table(index='SNP', columns='CLST', values='NCHROBS').fillna(0)
mac_minus_observed = mac_total - pivot
output = pivot.astype(int).astype(str) + ',' + mac_minus_observed.astype(int).astype(str)

with gzip.open('$OUTPUT_DIR/treemix/input.frq.gz', 'wt') as fh:
    fh.write(' '.join(pivot.columns) + '\n')
    for _, row in output.iterrows():
        fh.write(' '.join(row.values) + '\n')
PY

echo "[4/4] Run TreeMix with m=0..5 migration edges"
for m in 0 1 2 3 4 5; do
    if [ ! -f "$OUTPUT_DIR/treemix/output_m${m}.vertices.gz" ]; then
        treemix -i "$OUTPUT_DIR/treemix/input.frq.gz" \
            -m "$m" -bootstrap -k 1000 \
            -o "$OUTPUT_DIR/treemix/output_m${m}"
    fi
done

# OptM analysis (Evanno method). Use unquoted heredoc so bash expands $OUTPUT_DIR
# before R evaluates the string.
TREEMIX_DIR="$OUTPUT_DIR/treemix/"
Rscript - <<RS || true
if (require('OptM')) {
    optM_result <- optM(folder = '$TREEMIX_DIR', method = 'Evanno')
    print(optM_result)
} else {
    cat('Install OptM: install.packages("OptM")\n')
}
RS

echo "Done. Dsuite output: $OUTPUT_DIR/dsuite; TreeMix output: $OUTPUT_DIR/treemix"
