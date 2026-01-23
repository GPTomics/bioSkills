#!/bin/bash
# Metagene and heatmap visualization with deepTools

IP_BAM=$1
INPUT_BAM=$2
GENES_BED=$3
OUTPUT_PREFIX=${4:-"m6a_viz"}

# Create log2 IP/Input ratio bigWig
# Normalize by read depth
bamCompare \
    -b1 $IP_BAM \
    -b2 $INPUT_BAM \
    --scaleFactorsMethod readCount \
    --ratio log2 \
    --binSize 10 \
    -p 8 \
    -o ${OUTPUT_PREFIX}_log2ratio.bw

# Metagene around gene body
# Scale all genes to same length, add flanking regions
computeMatrix scale-regions \
    -S ${OUTPUT_PREFIX}_log2ratio.bw \
    -R $GENES_BED \
    --regionBodyLength 3000 \
    --beforeRegionStartLength 1000 \
    --afterRegionStartLength 1000 \
    --skipZeros \
    -p 8 \
    -o ${OUTPUT_PREFIX}_matrix.gz

# Profile plot (metagene)
plotProfile \
    -m ${OUTPUT_PREFIX}_matrix.gz \
    -o ${OUTPUT_PREFIX}_profile.pdf \
    --plotTitle "m6A signal across genes" \
    --yAxisLabel "log2(IP/Input)" \
    --regionsLabel "Genes"

# Heatmap
plotHeatmap \
    -m ${OUTPUT_PREFIX}_matrix.gz \
    -o ${OUTPUT_PREFIX}_heatmap.pdf \
    --colorMap RdBu_r \
    --whatToShow 'heatmap and colorbar' \
    --sortUsing mean

echo "Outputs: ${OUTPUT_PREFIX}_profile.pdf, ${OUTPUT_PREFIX}_heatmap.pdf"
