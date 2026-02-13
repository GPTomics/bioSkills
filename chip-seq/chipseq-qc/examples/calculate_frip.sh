#!/bin/bash
# Reference: MACS3 3.0+, Subread 2.0+, bedtools 2.31+, deepTools 3.5+, pybedtools 0.9+, pysam 0.22+, samtools 1.19+ | Verify API if version differs
# Calculate FRiP (Fraction of Reads in Peaks)

BAM=$1
PEAKS=$2

if [ -z "$BAM" ] || [ -z "$PEAKS" ]; then
    echo "Usage: $0 <bam_file> <peaks_file>"
    exit 1
fi

reads_in_peaks=$(bedtools intersect -a "$BAM" -b "$PEAKS" -u | samtools view -c -)
total_reads=$(samtools view -c -F 260 "$BAM")

frip=$(echo "scale=4; $reads_in_peaks / $total_reads" | bc)

echo "Total reads: $total_reads"
echo "Reads in peaks: $reads_in_peaks"
echo "FRiP: $frip"
