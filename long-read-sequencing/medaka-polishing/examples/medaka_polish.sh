#!/bin/bash

READS="reads.fastq.gz"
DRAFT="draft_assembly.fa"
OUTPUT_DIR="medaka_output"
THREADS=4
MODEL="r1041_e82_400bps_sup_v5.0.0"

medaka_consensus \
    -i $READS \
    -d $DRAFT \
    -o $OUTPUT_DIR \
    -t $THREADS \
    -m $MODEL

if [ -f "${OUTPUT_DIR}/consensus.fasta" ]; then
    echo "Polishing complete!"
    echo "Output: ${OUTPUT_DIR}/consensus.fasta"

    echo ""
    echo "Assembly statistics:"
    echo "Draft:"
    seqkit stats $DRAFT
    echo "Polished:"
    seqkit stats ${OUTPUT_DIR}/consensus.fasta
else
    echo "Error: Polishing failed"
    exit 1
fi
