#!/bin/bash
# Reference: minimap2 2.26+, samtools 1.19+, IsoQuant 3.5+, SQANTI3 5.2+, FLAIR 2.0+ | Verify CLI flags if version differs
# Long-read splicing analysis pipeline
#
# Workflow:
# 1. Splice-aware alignment with minimap2 (preset depends on platform)
# 2. Isoform discovery + quantification (IsoQuant for de novo, FLAIR for end-to-end)
# 3. Classification with SQANTI3 (FSM/ISM/NIC/NNC + artifact flags)
# 4. Differential analysis with FLAIR diffSplice or rMATS-long
#
# Set platform via PLATFORM=hifi or PLATFORM=ont

set -euo pipefail

PLATFORM=${PLATFORM:-hifi}
THREADS=${THREADS:-16}
REFERENCE=reference.fa
GTF=gencode.v45.annotation.gtf
SAMPLE=sample
FASTQ=${SAMPLE}.fastq.gz
OUTPUT_DIR=longread_output_${SAMPLE}

mkdir -p ${OUTPUT_DIR}

if [ "${PLATFORM}" = "hifi" ]; then
    PRESET="splice:hq"
    DATA_TYPE="pacbio_ccs"
elif [ "${PLATFORM}" = "ont" ]; then
    PRESET="splice -uf -k14"
    DATA_TYPE="nanopore"
else
    echo "PLATFORM must be hifi or ont" && exit 1
fi

# 1. Splice-aware alignment
minimap2 -ax ${PRESET} \
    -t ${THREADS} \
    --secondary=no \
    ${REFERENCE} ${FASTQ} | \
    samtools sort -@ ${THREADS} -o ${OUTPUT_DIR}/${SAMPLE}_aligned.bam
samtools index ${OUTPUT_DIR}/${SAMPLE}_aligned.bam

# 2a. Isoform discovery and quantification with IsoQuant (current SOTA)
isoquant.py \
    --reference ${REFERENCE} \
    --genedb ${GTF} \
    --bam ${OUTPUT_DIR}/${SAMPLE}_aligned.bam \
    --data_type ${DATA_TYPE} \
    --output ${OUTPUT_DIR}/isoquant \
    --threads ${THREADS} \
    --prefix ${SAMPLE} \
    --model_construction_strategy default_${DATA_TYPE}

# 2b. Alternative: FLAIR end-to-end pipeline
# Convert BAM to BED12 for FLAIR
bedtools bamtobed -bed12 -i ${OUTPUT_DIR}/${SAMPLE}_aligned.bam > ${OUTPUT_DIR}/${SAMPLE}.bed

flair correct \
    --query ${OUTPUT_DIR}/${SAMPLE}.bed \
    --genome ${REFERENCE} \
    --gtf ${GTF} \
    --output ${OUTPUT_DIR}/flair_corrected_${SAMPLE} \
    --threads ${THREADS}

flair collapse \
    --query ${OUTPUT_DIR}/flair_corrected_${SAMPLE}_all_corrected.bed \
    --reads ${FASTQ} \
    --genome ${REFERENCE} \
    --gtf ${GTF} \
    --output ${OUTPUT_DIR}/flair_collapsed_${SAMPLE} \
    --threads ${THREADS}

# 3. SQANTI3 classification (use SQANTI-LR branch for long reads)
sqanti3_qc.py \
    ${OUTPUT_DIR}/isoquant/${SAMPLE}/${SAMPLE}.transcript_models.gtf \
    ${GTF} \
    ${REFERENCE} \
    --output ${OUTPUT_DIR}/sqanti3 \
    --aligner_choice minimap2 \
    --cpus ${THREADS}

sqanti3_filter.py rules \
    ${OUTPUT_DIR}/sqanti3/${SAMPLE}_classification.txt \
    --gtf ${OUTPUT_DIR}/isoquant/${SAMPLE}/${SAMPLE}.transcript_models.gtf \
    --output ${OUTPUT_DIR}/sqanti3_filtered

echo "Pipeline complete. Outputs in ${OUTPUT_DIR}/"
echo "  IsoQuant transcript models: ${OUTPUT_DIR}/isoquant/${SAMPLE}/${SAMPLE}.transcript_models.gtf"
echo "  FLAIR isoforms: ${OUTPUT_DIR}/flair_collapsed_${SAMPLE}.isoforms.fa"
echo "  SQANTI3 classification: ${OUTPUT_DIR}/sqanti3/${SAMPLE}_classification.txt"
echo "  SQANTI3 filtered (rt-switching, intra-priming removed): ${OUTPUT_DIR}/sqanti3_filtered"
