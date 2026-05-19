#!/usr/bin/env bash
# Reference: TOGA 1.1.7+, CESAR 2.0+, Cactus 2.9+, HAL toolkit 2.3+, Nextflow 24+, UCSC kentUtils 2024+ | Verify API if version differs
# TOGA annotation projection from a Cactus HAL alignment.

set -euo pipefail

CACTUS_HAL=${1:?usage: $0 CACTUS_HAL REFERENCE_NAME QUERY_NAME REF_ANNOTATION_BED OUTPUT_DIR}
REFERENCE_NAME=${2:?missing reference genome name (e.g. GRCh38)}
QUERY_NAME=${3:?missing query genome name}
REF_ANNOTATION_BED=${4:?missing reference annotation BED12 file}
OUTPUT_DIR=${5:?missing output dir}
THREADS=${THREADS:-32}

mkdir -p "$OUTPUT_DIR"/{chain,toga}

echo "[1/4] Extract syntenic blocks from HAL"
halSynteny "$CACTUS_HAL" "$REFERENCE_NAME" "$QUERY_NAME" \
    > "$OUTPUT_DIR/chain/syntenic_blocks.bed"

echo "[2/4] Convert HAL to UCSC chain format"
# Note: TOGA expects UCSC-style chain files
# halSynteny -> bed -> chain (manual conversion or via nf-core/toga workflow)
hal2maf "$CACTUS_HAL" "$REFERENCE_NAME" \
    --refGenome "$REFERENCE_NAME" \
    --chunkSize 1000000 \
    > "$OUTPUT_DIR/chain/${QUERY_NAME}.maf"

# Convert MAF to chain
# This requires UCSC kentUtils mafChain or similar wrapper
# In production, use the TOGA Nextflow pipeline which handles this

echo "[3/4] Prepare 2bit files for TOGA"
# TOGA expects .2bit reference + query genome files
faToTwoBit "$REFERENCE_NAME.fa" "$OUTPUT_DIR/chain/${REFERENCE_NAME}.2bit"
faToTwoBit "$QUERY_NAME.fa" "$OUTPUT_DIR/chain/${QUERY_NAME}.2bit"

echo "[4/4] Run TOGA Nextflow pipeline"
nextflow run hillerlab/TOGA \
    --chain "$OUTPUT_DIR/chain/${QUERY_NAME}.chain.gz" \
    --bed "$REF_ANNOTATION_BED" \
    --tDB "$OUTPUT_DIR/chain/${REFERENCE_NAME}.2bit" \
    --qDB "$OUTPUT_DIR/chain/${QUERY_NAME}.2bit" \
    --nextflow_dir "$OUTPUT_DIR/toga/nextflow_run" \
    --pn "${QUERY_NAME}_TOGA" \
    --cpu "$THREADS" \
    --quiet \
    --output "$OUTPUT_DIR/toga"

# Output:
#   $OUTPUT_DIR/toga/${QUERY_NAME}_TOGA/loss_summ_data.tsv      Per-gene intactness
#   $OUTPUT_DIR/toga/${QUERY_NAME}_TOGA/orthology_classification.tsv  Orthology calls
#   $OUTPUT_DIR/toga/${QUERY_NAME}_TOGA/query_annotation.bed    Projected annotation
#   $OUTPUT_DIR/toga/${QUERY_NAME}_TOGA/query_annotation.gff    GFF format

echo "Parse TOGA intactness summary"
python3 - <<PY
import pandas as pd
import os

summary = '$OUTPUT_DIR/toga/${QUERY_NAME}_TOGA/loss_summ_data.tsv'
if os.path.exists(summary):
    df = pd.read_csv(summary, sep='\t')
    # Standard columns: TRANSCRIPT, STATUS, IS_INTACT, ...
    if 'STATUS' in df.columns:
        counts = df['STATUS'].value_counts()
        print('TOGA intactness distribution (I/PI/UL/L/M/PM):')
        print(counts)
        intact = df[df['STATUS'] == 'I']
        partial = df[df['STATUS'] == 'PI']
        print(f'\nIntact (I): {len(intact)} / {len(df)} ({len(intact)/len(df)*100:.1f}%)')
        print(f'Partial Intact (PI): {len(partial)} / {len(df)} ({len(partial)/len(df)*100:.1f}%)')
PY

echo "Done. TOGA output: $OUTPUT_DIR/toga/${QUERY_NAME}_TOGA/"
