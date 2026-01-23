#!/bin/bash
# Iso-Seq3 pipeline for isoform discovery

SUBREADS=$1
PRIMERS=$2
OUTPUT_DIR=${3:-"isoseq_output"}
THREADS=${4:-8}

mkdir -p $OUTPUT_DIR

# Step 1: Generate CCS reads
# min-rq 0.9: Minimum read quality (standard for Iso-Seq)
# min-passes 3: Minimum polymerase passes
echo "Generating CCS reads..."
ccs $SUBREADS ${OUTPUT_DIR}/ccs.bam \
    --min-rq 0.9 \
    --min-passes 3 \
    --num-threads $THREADS

# Step 2: Primer removal
# --isoseq: Iso-Seq specific settings
echo "Removing primers..."
lima ${OUTPUT_DIR}/ccs.bam $PRIMERS ${OUTPUT_DIR}/demux.bam \
    --isoseq \
    --num-threads $THREADS

# Find demux output (depends on primer names)
DEMUX=$(ls ${OUTPUT_DIR}/demux.*--*.bam | head -1)

# Step 3: Refine
# Removes polyA tails and concatemers
echo "Refining reads..."
isoseq3 refine $DEMUX $PRIMERS ${OUTPUT_DIR}/refined.bam \
    --require-polya

# Step 4: Cluster into isoforms
# Generates high-quality (HQ) and low-quality (LQ) isoforms
echo "Clustering isoforms..."
isoseq3 cluster ${OUTPUT_DIR}/refined.bam ${OUTPUT_DIR}/clustered.bam \
    --verbose \
    --num-threads $THREADS

echo "Pipeline complete."
echo "HQ transcripts: ${OUTPUT_DIR}/clustered.hq.fasta.gz"
echo "LQ transcripts: ${OUTPUT_DIR}/clustered.lq.fasta.gz"
echo "Cluster report: ${OUTPUT_DIR}/clustered.cluster_report.csv"
