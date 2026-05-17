#!/bin/bash
# Reference: deepTools 3.5+, samtools 1.19+ | Verify API if version differs
# deepTools heatmap pipeline: bamCoverage CPM-normalized bigWig -> TSS-centered
# matrix -> heatmap + average profile. Usage: edit BAM_FILE/GENES_BED below.
set -euo pipefail

BAM_FILE="${1:-sample.bam}"
OUTPUT_DIR="${2:-visualization}"
GENES_BED="${3:-genes.bed}"

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
