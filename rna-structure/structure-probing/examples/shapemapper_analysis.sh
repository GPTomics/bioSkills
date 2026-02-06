#!/bin/bash
# SHAPE-MaP analysis with ShapeMapper2.
# Processes modified and untreated samples to generate per-nucleotide reactivities.

TARGET_FA=$1
NAME=${2:-"my_rna"}
MOD_R1=$3
MOD_R2=$4
UNMOD_R1=$5
UNMOD_R2=$6
OUTPUT_DIR=${7:-"shapemapper_results"}
THREADS=${8:-8}

if [ -z "$TARGET_FA" ] || [ -z "$MOD_R1" ]; then
    echo "Usage: $0 <target.fa> <name> <mod_R1.fq.gz> <mod_R2.fq.gz> <unmod_R1.fq.gz> <unmod_R2.fq.gz> [output_dir] [threads]"
    echo ""
    echo "Runs ShapeMapper2 on SHAPE-MaP data to generate reactivity profiles."
    echo "Note: ShapeMapper2 is Linux-only. Use Docker on macOS:"
    echo "  docker run -v \$(pwd):/data shapemapper2/shapemapper2 shapemapper \\"
    echo "    --target /data/target.fa --modified --R1 /data/mod_R1.fq.gz ..."
    exit 1
fi

echo "=== Running ShapeMapper2 ==="
echo "Target: $TARGET_FA"
echo "Name: $NAME"
echo "Output: $OUTPUT_DIR"

shapemapper \
    --target "$TARGET_FA" \
    --name "$NAME" \
    --modified --R1 "$MOD_R1" --R2 "$MOD_R2" \
    --untreated --R1 "$UNMOD_R1" --R2 "$UNMOD_R2" \
    --out "$OUTPUT_DIR" \
    --nproc "$THREADS" \
    --min-depth 5000 \
    --overwrite

echo ""
echo "=== QC Summary ==="
if [ -f "$OUTPUT_DIR/${NAME}_shapemapper_log.txt" ]; then
    grep -E "Effective read depth|Mutation rate|quality" "$OUTPUT_DIR/${NAME}_shapemapper_log.txt"
fi

echo ""
echo "=== Output files ==="
echo "  ${OUTPUT_DIR}/${NAME}_profile.txt        - Reactivity profile"
echo "  ${OUTPUT_DIR}/${NAME}_map.shape           - SHAPE file for RNAfold"
echo "  ${OUTPUT_DIR}/${NAME}_profile.pdf         - Reactivity bar plot"

echo ""
echo "=== Use reactivities to constrain folding ==="
echo "  RNAfold --shape=${OUTPUT_DIR}/${NAME}_map.shape < ${TARGET_FA}"
