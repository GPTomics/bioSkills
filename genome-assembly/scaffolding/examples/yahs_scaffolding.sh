#!/bin/bash
set -euo pipefail

ASSEMBLY=$1
HIC_R1=$2
HIC_R2=$3
OUTPUT_PREFIX=${4:-scaffolds}
THREADS=${5:-16}

echo "Indexing assembly..."
samtools faidx "$ASSEMBLY"
bwa index "$ASSEMBLY"

echo "Aligning Hi-C reads..."
bwa mem -5SP -t "$THREADS" "$ASSEMBLY" "$HIC_R1" "$HIC_R2" | \
    samtools view -@ 8 -bhS - | \
    samtools sort -@ 8 -n -o aligned_sorted.bam

echo "Converting to BED format..."
bedtools bamtobed -i aligned_sorted.bam | \
    sort -k4 > aligned.bed

echo "Running YaHS scaffolding..."
yahs "$ASSEMBLY" aligned.bed -o "$OUTPUT_PREFIX"

SCAFFOLDS="${OUTPUT_PREFIX}_scaffolds_final.fa"

echo "Generating assembly statistics..."
seqkit stats "$SCAFFOLDS" > "${OUTPUT_PREFIX}_stats.txt"
echo "Original assembly:"
seqkit stats "$ASSEMBLY" >> "${OUTPUT_PREFIX}_stats.txt"

echo "Creating contact map for visualization..."
samtools faidx "$SCAFFOLDS"

juicer pre "${OUTPUT_PREFIX}.bin" "${OUTPUT_PREFIX}_scaffolds_final.agp" \
    "${ASSEMBLY}.fai" 2>/dev/null | \
    sort -k2,2d -k6,6d -T ./ --parallel="$THREADS" -S 32G | \
    awk 'NF' > "${OUTPUT_PREFIX}.pre.txt"

java -Xmx32G -jar juicer_tools.jar pre \
    "${OUTPUT_PREFIX}.pre.txt" \
    "${OUTPUT_PREFIX}.hic" \
    <(cut -f1,2 "${SCAFFOLDS}.fai")

echo "Running BUSCO..."
busco -i "$SCAFFOLDS" -l eukaryota_odb10 -o "${OUTPUT_PREFIX}_busco" \
    -m genome -c "$THREADS" --quiet

echo "Done! Output files:"
echo "  Scaffolds: $SCAFFOLDS"
echo "  Contact map: ${OUTPUT_PREFIX}.hic"
echo "  Statistics: ${OUTPUT_PREFIX}_stats.txt"
