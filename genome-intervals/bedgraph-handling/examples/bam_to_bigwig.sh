#!/bin/bash
# Convert BAM to normalized bigWig via bedGraph

BAM=$1
CHROM_SIZES=$2

if [ -z "$BAM" ] || [ -z "$CHROM_SIZES" ]; then
    echo "Usage: $0 <bam_file> <chrom.sizes>"
    exit 1
fi

NAME=$(basename $BAM .bam)

total_reads=$(samtools view -c -F 260 $BAM)
scale=$(echo "scale=10; 1000000 / $total_reads" | bc)
echo "Total reads: $total_reads, Scale factor: $scale"

echo "Generating bedGraph..."
bedtools genomecov -ibam $BAM -bg -scale $scale > ${NAME}.tmp.bedgraph

echo "Sorting..."
sort -k1,1 -k2,2n ${NAME}.tmp.bedgraph > ${NAME}.sorted.bedgraph

echo "Clipping to chromosome boundaries..."
bedClip ${NAME}.sorted.bedgraph $CHROM_SIZES ${NAME}.clipped.bedgraph

echo "Converting to bigWig..."
bedGraphToBigWig ${NAME}.clipped.bedgraph $CHROM_SIZES ${NAME}.bw

rm -f ${NAME}.tmp.bedgraph ${NAME}.sorted.bedgraph ${NAME}.clipped.bedgraph

echo "Done: ${NAME}.bw"
