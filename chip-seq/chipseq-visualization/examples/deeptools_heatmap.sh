#!/bin/bash

BAM_FILE="sample.bam"
OUTPUT_DIR="visualization"
GENES_BED="genes.bed"

mkdir -p $OUTPUT_DIR

bamCoverage \
    -b $BAM_FILE \
    -o ${OUTPUT_DIR}/sample.bw \
    --normalizeUsing CPM \
    --binSize 10 \
    --numberOfProcessors 8

computeMatrix reference-point \
    --referencePoint TSS \
    -b 3000 -a 3000 \
    -R $GENES_BED \
    -S ${OUTPUT_DIR}/sample.bw \
    -o ${OUTPUT_DIR}/matrix_tss.gz \
    --numberOfProcessors 8

plotHeatmap \
    -m ${OUTPUT_DIR}/matrix_tss.gz \
    -o ${OUTPUT_DIR}/heatmap_tss.png \
    --colorMap RdBu_r \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --refPointLabel TSS \
    --plotTitle 'ChIP-seq Signal at TSS'

plotProfile \
    -m ${OUTPUT_DIR}/matrix_tss.gz \
    -o ${OUTPUT_DIR}/profile_tss.png \
    --plotTitle 'Average Signal Profile' \
    --refPointLabel TSS

echo "Visualization complete. Output in ${OUTPUT_DIR}/"
