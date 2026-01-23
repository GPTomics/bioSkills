#!/bin/bash
# m6A peak calling with MACS3

IP_BAMS="IP_rep1.bam IP_rep2.bam"
INPUT_BAMS="Input_rep1.bam Input_rep2.bam"
OUTPUT_PREFIX="m6a_peaks"
GENOME_SIZE="hs"  # hs for human, mm for mouse

# Call peaks
# --nomodel: Don't build shifting model (m6A peaks differ from TF ChIP)
# --extsize 150: Extend reads to ~150bp (typical fragment size)
# -q 0.05: FDR threshold
macs3 callpeak \
    -t $IP_BAMS \
    -c $INPUT_BAMS \
    -f BAMPE \
    -g $GENOME_SIZE \
    -n $OUTPUT_PREFIX \
    --nomodel \
    --extsize 150 \
    -q 0.05 \
    --keep-dup auto

# Filter peaks by fold enrichment
# FC > 2: Require 2-fold IP/Input enrichment
awk -F'\t' 'NR > 1 && $7 > 2' ${OUTPUT_PREFIX}_peaks.xls > ${OUTPUT_PREFIX}_filtered.bed

echo "Peaks called: $(wc -l < ${OUTPUT_PREFIX}_peaks.narrowPeak)"
echo "After FC>2 filter: $(wc -l < ${OUTPUT_PREFIX}_filtered.bed)"
