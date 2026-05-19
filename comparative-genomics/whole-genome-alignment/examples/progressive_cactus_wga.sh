#!/usr/bin/env bash
# Reference: Cactus 2.9+, HAL toolkit 2.3+, RepeatMasker 4.1.5+, RepeatModeler 2.0.5+, Toil 6.0+ | Verify API if version differs
# Progressive Cactus reference-free multi-genome alignment pipeline.
# Demonstrates: per-species softmasking, seqFile preparation, Cactus invocation, HAL inspection.

set -euo pipefail

GENOMES_DIR=${1:?usage: $0 GENOMES_DIR SPECIES_TREE OUTPUT_DIR}
SPECIES_TREE=${2:?missing species tree (Newick with branch lengths)}
OUTPUT_DIR=${3:?missing output dir}
THREADS=${THREADS:-16}
BRANCH_SCALE=${BRANCH_SCALE:-1.0}  # 1.0 vertebrate; 0.5-0.7 plant; 0.3-0.5 bacteria

mkdir -p "$OUTPUT_DIR"/{masked,cactus_run}

echo "[1/4] Mask repeats per species"
for fa in "$GENOMES_DIR"/*.fa; do
    species=$(basename "$fa" .fa)
    if [ ! -f "$OUTPUT_DIR/masked/${species}.fa.masked" ]; then
        BuildDatabase -name "$OUTPUT_DIR/masked/${species}_DB" "$fa"
        RepeatModeler -database "$OUTPUT_DIR/masked/${species}_DB" -threads "$THREADS"
        RepeatMasker -lib "$OUTPUT_DIR/masked/${species}_DB-families.fa" \
            -xsmall -pa "$THREADS" -dir "$OUTPUT_DIR/masked" "$fa"
    fi
done

echo "[2/4] Build Cactus seqFile (tree + paths)"
SEQFILE="$OUTPUT_DIR/seqFile.txt"
cat "$SPECIES_TREE" > "$SEQFILE"
echo "" >> "$SEQFILE"
for fa in "$OUTPUT_DIR/masked"/*.fa.masked; do
    species=$(basename "$fa" .fa.masked)
    echo -e "${species}\t$(realpath "$fa")" >> "$SEQFILE"
done

echo "[3/4] Run Progressive Cactus"
JOBSTORE="$OUTPUT_DIR/cactus_run/jobStore"
HAL_OUT="$OUTPUT_DIR/cactus_run/output.hal"

if [ ! -f "$HAL_OUT" ]; then
    cactus "$JOBSTORE" "$SEQFILE" "$HAL_OUT" \
        --binariesMode local \
        --workDir "$OUTPUT_DIR/cactus_run/work" \
        --maxCores "$THREADS" \
        --branchScale "$BRANCH_SCALE" \
        --logFile "$OUTPUT_DIR/cactus_run/cactus.log"
fi

echo "[4/4] Inspect HAL output"
halStats "$HAL_OUT"

# Extract syntenic blocks for downstream TOGA / annotation projection
if command -v halSynteny &> /dev/null; then
    REF_GENOME=$(awk 'NR>1{print $1; exit}' "$SEQFILE")
    halSynteny "$HAL_OUT" "$REF_GENOME" \
        --refSeq chr1 > "$OUTPUT_DIR/cactus_run/syntenic_blocks_chr1.bed"
fi

echo "Done. HAL: $HAL_OUT"
