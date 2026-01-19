#!/bin/bash
# Basecall with Dorado

INPUT_DIR=$1
OUTPUT_BAM=${2:-calls.bam}
MODEL=${3:-sup}

if [ -z "$INPUT_DIR" ]; then
    echo "Usage: $0 <pod5_dir> [output.bam] [model]"
    echo "Models: fast, hac, sup (default: sup)"
    exit 1
fi

echo "Basecalling $INPUT_DIR with $MODEL model..."
dorado basecaller $MODEL $INPUT_DIR > $OUTPUT_BAM

echo "Generating FASTQ..."
samtools fastq $OUTPUT_BAM | gzip > ${OUTPUT_BAM%.bam}.fastq.gz

echo "Stats:"
samtools stats $OUTPUT_BAM | grep "^SN" | head -10

echo "Done: $OUTPUT_BAM"
