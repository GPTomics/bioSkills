#!/bin/bash
# Comprehensive alignment validation

BAM=$1

if [ -z "$BAM" ]; then
    echo "Usage: $0 <bam_file>"
    exit 1
fi

NAME=$(basename $BAM .bam)

echo "=== Alignment Validation: $NAME ==="

echo -e "\n--- Basic Stats ---"
samtools flagstat $BAM

echo -e "\n--- Mapping Rate ---"
mapped=$(samtools view -c -F 4 $BAM)
total=$(samtools view -c $BAM)
echo "Mapped: $mapped / $total ($(echo "scale=1; $mapped/$total*100" | bc)%)"

echo -e "\n--- Proper Pairing ---"
proper=$(samtools view -c -f 2 $BAM)
echo "Properly paired: $proper ($(echo "scale=1; $proper/$mapped*100" | bc)%)"

echo -e "\n--- Insert Size ---"
samtools stats $BAM 2>/dev/null | grep "insert size" | head -3

echo -e "\n--- Strand Balance ---"
fwd=$(samtools view -c -F 16 $BAM)
rev=$(samtools view -c -f 16 $BAM)
echo "Forward: $fwd, Reverse: $rev, Ratio: $(echo "scale=3; $fwd/$rev" | bc)"

echo -e "\n--- MAPQ Distribution ---"
echo "MAPQ 0 (multi-mapper): $(samtools view -q 0 -Q 1 $BAM | wc -l)"
echo "MAPQ >= 30: $(samtools view -c -q 30 $BAM)"

echo -e "\n--- Top Chromosomes ---"
samtools idxstats $BAM | sort -k3 -rn | head -10
