#!/usr/bin/env bash
# Reference: deepTools 3.5+, samtools 1.19+ | Verify with `deeptools --version`, `samtools --version` if installed releases differ.
# Genome-coordinate metagene + peak-centred heatmap via deepTools. Pairs with Guitar transcript-feature metagene.
# deepTools `--operation log2` is the modern syntax; `--ratio log2` is being phased out.

set -euo pipefail

IP_BAM=${1:-'aligned/IP_rep1_Aligned.sortedByCoord.out.bam'}
INPUT_BAM=${2:-'aligned/Input_rep1_Aligned.sortedByCoord.out.bam'}
GENES_BED=${3:-'refs/protein_coding.bed'}
PEAKS_BED=${4:-'exomepeak2_output/m6a_run1/peaks.bed'}
OUTDIR=${5:-'figures'}
THREADS=8

mkdir -p "${OUTDIR}"

# Step 1: IP-over-Input log2 bigWig.
bamCompare \
    -b1 "${IP_BAM}" \
    -b2 "${INPUT_BAM}" \
    --operation log2 \
    --pseudocount 1 \
    --binSize 25 \
    --normalizeUsing CPM \
    --numberOfProcessors "${THREADS}" \
    -o "${OUTDIR}/log2_IP_over_Input.bw"

# Step 2: Genome-coordinate metagene over protein-coding genes (5' end to 3' end, scaled).
computeMatrix scale-regions \
    --regionsFileName "${GENES_BED}" \
    --scoreFileName "${OUTDIR}/log2_IP_over_Input.bw" \
    --regionBodyLength 2000 \
    --upstream 500 \
    --downstream 500 \
    --skipZeros \
    --numberOfProcessors "${THREADS}" \
    --outFileName "${OUTDIR}/genebody_matrix.gz"

plotProfile \
    --matrixFile "${OUTDIR}/genebody_matrix.gz" \
    --outFileName "${OUTDIR}/genebody_profile.pdf" \
    --plotTitle 'm6A IP/Input metagene (genome-coordinate, scaled gene body)' \
    --plotType lines

plotHeatmap \
    --matrixFile "${OUTDIR}/genebody_matrix.gz" \
    --outFileName "${OUTDIR}/genebody_heatmap.pdf" \
    --colorMap RdBu_r \
    --plotTitle 'm6A IP/Input across gene body'

# Step 3: Peak-centred matrix + heatmap (+/-500 bp around peak centre).
computeMatrix reference-point \
    --regionsFileName "${PEAKS_BED}" \
    --scoreFileName "${OUTDIR}/log2_IP_over_Input.bw" \
    --referencePoint center \
    --upstream 500 \
    --downstream 500 \
    --binSize 25 \
    --skipZeros \
    --numberOfProcessors "${THREADS}" \
    --outFileName "${OUTDIR}/peak_centred_matrix.gz" \
    --outFileNameMatrix "${OUTDIR}/peak_centred_matrix.tab"

plotHeatmap \
    --matrixFile "${OUTDIR}/peak_centred_matrix.gz" \
    --outFileName "${OUTDIR}/peak_centred_heatmap.pdf" \
    --kmeans 3 \
    --colorMap viridis \
    --plotTitle 'm6A signal centred at peaks (k-means k=3)'

echo "Outputs:"
echo "  - ${OUTDIR}/genebody_profile.pdf (genome-coordinate metagene profile)"
echo "  - ${OUTDIR}/genebody_heatmap.pdf (genome-coordinate heatmap)"
echo "  - ${OUTDIR}/peak_centred_heatmap.pdf (peak-centred k-means k=3)"
echo "Pair with Guitar transcript-feature metagene for 5UTR/CDS/3UTR semantics."
