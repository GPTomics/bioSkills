#!/usr/bin/env bash
# Reference: skani 0.2.5+, GTDB-Tk 2.7+, CheckM2 1.0.2+, FastANI 1.34+, Python 3.10+, pandas 2.2+ | Verify API if version differs
# skani + GTDB-Tk pipeline for genome-distance and species delineation.

set -euo pipefail

GENOMES_DIR=${1:?usage: $0 GENOMES_DIR OUTPUT_DIR}
OUTPUT_DIR=${2:?missing output dir}
THREADS=${THREADS:-16}
GTDBTK_DATA_PATH=${GTDBTK_DATA_PATH:?must set GTDBTK_DATA_PATH to GTDB release directory}

mkdir -p "$OUTPUT_DIR"/{checkm2,skani,gtdbtk}

echo "[1/4] CheckM2 quality filtering (MAGs)"
checkm2 predict --threads "$THREADS" \
    --input "$GENOMES_DIR" \
    --output-directory "$OUTPUT_DIR/checkm2"

# Filter to >=70% completeness, <5% contamination
awk 'NR > 1 && $2 >= 70 && $3 < 5 {print $1}' \
    "$OUTPUT_DIR/checkm2/quality_report.tsv" > "$OUTPUT_DIR/checkm2/passed_genomes.txt"

echo "Passed: $(wc -l < "$OUTPUT_DIR/checkm2/passed_genomes.txt") / $(ls "$GENOMES_DIR"/*.fa | wc -l)"

echo "[2/4] skani all-vs-all ANI"
mkdir -p "$OUTPUT_DIR/skani/sketches"
# Sketch all passed genomes
while read name; do
    fa_path=$(find "$GENOMES_DIR" -name "${name}*" | head -1)
    if [ -n "$fa_path" ]; then
        cp "$fa_path" "$OUTPUT_DIR/skani/sketches/"
    fi
done < "$OUTPUT_DIR/checkm2/passed_genomes.txt"

skani triangle "$OUTPUT_DIR/skani/sketches"/*.fa \
    -t "$THREADS" --robust --sparse -o "$OUTPUT_DIR/skani/ani_matrix.tsv"
# --sparse emits Ref_file<TAB>Query_file<TAB>ANI<TAB>Align_fraction_ref<TAB>Align_fraction_query
# (the Python parser below depends on this columnar format).

echo "[3/4] Apply 95% ANI + AF >=0.5 species delineation"
python3 - <<PY
import pandas as pd

df = pd.read_csv('$OUTPUT_DIR/skani/ani_matrix.tsv', sep='\t')
# skani triangle --sparse columns: Ref_file Query_file ANI Align_fraction_ref Align_fraction_query
# Take the first 5 columns and rename for consistency; ignore any extra metadata columns.
df = df.iloc[:, :5].copy()
df.columns = ['ref', 'query', 'ani', 'af_ref', 'af_query']
df['min_af'] = df[['af_ref', 'af_query']].min(axis=1)

# Standard delineation
same_species = df[(df['ani'] >= 95.0) & (df['min_af'] >= 0.5)]
print(f'\n95% ANI + 50% AF same-species pairs: {len(same_species)}')

# Borderline cases for cross-validation
borderline = df[(df['ani'] >= 93) & (df['ani'] < 96) & (df['min_af'] >= 0.5)]
print(f'\nBorderline (93-96% ANI): {len(borderline)} pairs (cross-validate with OrthoANI)')

# Output cluster summary
same_species.to_csv('$OUTPUT_DIR/skani/same_species_pairs.tsv', sep='\t', index=False)
borderline.to_csv('$OUTPUT_DIR/skani/borderline_pairs.tsv', sep='\t', index=False)
PY

echo "[4/4] GTDB-Tk taxonomic classification"
gtdbtk classify_wf \
    --genome_dir "$OUTPUT_DIR/skani/sketches" \
    --out_dir "$OUTPUT_DIR/gtdbtk" \
    --cpus "$THREADS" \
    --extension fa \
    --skip_ani_screen  # optional; skip if want phylogeny-based only

# Parse GTDB taxonomy
if [ -f "$OUTPUT_DIR/gtdbtk/classify/gtdbtk.bac120.summary.tsv" ]; then
    python3 - <<PY
import pandas as pd

df = pd.read_csv('$OUTPUT_DIR/gtdbtk/classify/gtdbtk.bac120.summary.tsv', sep='\t')
classified = df['classification'].notna().sum()
species_level = df['classification'].str.contains('s__').sum() if 'classification' in df.columns else 0
print(f'\nGTDB classification: {classified} bacterial genomes')
print(f'Species-level: {species_level}')
print(df[['user_genome', 'classification', 'fastani_ani']].head(10))
PY
fi

echo "Done. ANI: $OUTPUT_DIR/skani; Taxonomy: $OUTPUT_DIR/gtdbtk"
