#!/bin/bash
# Complete basecalling pipeline

INPUT=$1
OUTPUT=${2:-basecall_output}
MODEL=${3:-hac}
MIN_QUAL=${4:-10}
MIN_LEN=${5:-500}

if [ -z "$INPUT" ]; then
    echo "Usage: $0 <input_dir> [output_dir] [model] [min_qual] [min_len]"
    exit 1
fi

mkdir -p $OUTPUT

if ls $INPUT/*.fast5 1> /dev/null 2>&1; then
    echo "Converting FAST5 to POD5..."
    mkdir -p $OUTPUT/pod5
    pod5 convert fast5 $INPUT/*.fast5 --output $OUTPUT/pod5/
    POD5_DIR="$OUTPUT/pod5"
elif ls $INPUT/*.pod5 1> /dev/null 2>&1; then
    POD5_DIR="$INPUT"
else
    echo "No FAST5 or POD5 files found in $INPUT"
    exit 1
fi

echo "Basecalling with $MODEL model..."
dorado basecaller $MODEL $POD5_DIR > $OUTPUT/calls.bam

echo "Converting to FASTQ..."
samtools fastq $OUTPUT/calls.bam | gzip > $OUTPUT/calls.fastq.gz

echo "Filtering (Q>=$MIN_QUAL, len>=$MIN_LEN)..."
gunzip -c $OUTPUT/calls.fastq.gz | \
    chopper -q $MIN_QUAL -l $MIN_LEN | \
    gzip > $OUTPUT/filtered.fastq.gz

echo "Generating QC report..."
NanoPlot --fastq $OUTPUT/filtered.fastq.gz -o $OUTPUT/qc/ --plots hex dot

echo "=== Summary ==="
echo "Raw reads: $(samtools view -c $OUTPUT/calls.bam)"
echo "Filtered reads: $(zcat $OUTPUT/filtered.fastq.gz | echo $((`wc -l`/4)))"
echo "QC report: $OUTPUT/qc/NanoPlot-report.html"
echo "Done!"
