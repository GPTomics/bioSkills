#!/bin/bash
# Reference: ChromHMM 1.27+, Java 8+, samtools 1.19+ | Verify API if version differs
# ChromHMM end-to-end pipeline: cellMarkFileTable -> BinarizeBam -> LearnModel
# at multiple N -> compare emission heatmaps -> OverlapEnrichment +
# NeighborhoodEnrichment. Canonical 15-state Roadmap-compatible workflow.

set -euo pipefail

CELL_TYPE=${1:-GM12878}
CHROMHMM_JAR=${2:-ChromHMM.jar}
CHROMSIZES=${3:-hg38.chromsizes.txt}
BAM_DIR=${4:-bams/}
ANCHOR_DIR=${5:-anchors/}     # BED files for OverlapEnrichment (CGI, repeats, etc.)
TSS_FILE=${6:-tss_hg38.txt}   # for NeighborhoodEnrichment


# === 1. Build cellMarkFileTable ===
# Format: cell_type<TAB>mark<TAB>file<TAB>control_file
# Use sonicated input as control (NOT IgG) for histone marks.
cat > cellMarkFileTable.txt << EOF
${CELL_TYPE}	H3K4me3	${BAM_DIR}/${CELL_TYPE}_h3k4me3.bam	${BAM_DIR}/${CELL_TYPE}_input.bam
${CELL_TYPE}	H3K27ac	${BAM_DIR}/${CELL_TYPE}_h3k27ac.bam	${BAM_DIR}/${CELL_TYPE}_input.bam
${CELL_TYPE}	H3K4me1	${BAM_DIR}/${CELL_TYPE}_h3k4me1.bam	${BAM_DIR}/${CELL_TYPE}_input.bam
${CELL_TYPE}	H3K36me3	${BAM_DIR}/${CELL_TYPE}_h3k36me3.bam	${BAM_DIR}/${CELL_TYPE}_input.bam
${CELL_TYPE}	H3K27me3	${BAM_DIR}/${CELL_TYPE}_h3k27me3.bam	${BAM_DIR}/${CELL_TYPE}_input.bam
${CELL_TYPE}	H3K9me3	${BAM_DIR}/${CELL_TYPE}_h3k9me3.bam	${BAM_DIR}/${CELL_TYPE}_input.bam
EOF


# === 2. Binarize BAMs ===
# Default bin size 200 bp; reduce to 100 or 50 only for sharp-boundary studies.
# Output: per-chromosome _binary.txt files in binarized_${CELL_TYPE}/
mkdir -p binarized_${CELL_TYPE}
java -mx16G -jar $CHROMHMM_JAR BinarizeBam \
    -b 200 \
    -peaks \
    $CHROMSIZES \
    $BAM_DIR \
    cellMarkFileTable.txt \
    binarized_${CELL_TYPE}/


# === 3. Learn model at multiple N — compare emission matrices ===
# Train 15-, 18-, and 25-state models; choose smallest N where states are
# biologically distinct (inspect emissions_N.png).
for N in 15 18 25; do
    java -mx16G -jar $CHROMHMM_JAR LearnModel \
        -p 8 \
        -printposterior \
        binarized_${CELL_TYPE}/ \
        model_${N}state_${CELL_TYPE}/ \
        $N \
        hg38
    # Output: model_${N}.txt (emission + transition matrices)
    #         emissions_${N}.png, transitions_${N}.png
    #         ${CELL_TYPE}_${N}_segments.bed (state assignments)
done


# === 4. OverlapEnrichment with anchor features ===
# Anchor BED files: CpG islands, RefSeq exons, lncRNA, repeats, blacklist, etc.
# Typically downloaded from UCSC tracks. ChromHMM iterates over all .bed in dir.
java -mx16G -jar $CHROMHMM_JAR OverlapEnrichment \
    -labels \
    model_15state_${CELL_TYPE}/${CELL_TYPE}_15_segments.bed \
    $ANCHOR_DIR \
    overlap_enrichment_15state_${CELL_TYPE}
# Output: per-state enrichment over each anchor BED


# === 5. NeighborhoodEnrichment around TSS ===
# Anchor TSS positions (one position per gene). Reveals state distribution
# relative to TSS — useful for verifying state interpretations.
java -mx16G -jar $CHROMHMM_JAR NeighborhoodEnrichment \
    -labels \
    model_15state_${CELL_TYPE}/${CELL_TYPE}_15_segments.bed \
    $TSS_FILE \
    neighborhood_TSS_15state_${CELL_TYPE}


# === 6. Apply Roadmap precomputed 25-state model (optional) ===
# Use for cross-cell-type compatibility with Roadmap Epigenomics annotations.
# Note: mark panel must match Roadmap's exactly (5 core marks).
# wget http://compbio.mit.edu/ChromHMM/model_25_imputed12marks.txt
# java -mx16G -jar $CHROMHMM_JAR MakeSegmentation \
#     model_25_imputed12marks.txt \
#     binarized_${CELL_TYPE}/ \
#     roadmap_25state_segmentation_${CELL_TYPE}


# === 7. Summary ===
echo "=== Summary: $CELL_TYPE ==="
for N in 15 18 25; do
    BED=model_${N}state_${CELL_TYPE}/${CELL_TYPE}_${N}_segments.bed
    [ -f $BED ] || continue
    echo "$N states: $(cut -f4 $BED | sort -u | wc -l) unique state labels"
    echo "  Bins per state:"
    cut -f4 $BED | sort | uniq -c | sort -nr | head -10
done
