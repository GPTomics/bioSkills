#!/bin/bash
# Reference: Bismark 0.24+, Bowtie2 2.5+, Trim Galore 0.6.10+, samtools 1.19+ | Verify API if version differs
#
# End-to-end WGBS/EM-seq paired-end alignment: prepare index -> trim -> align ->
# deduplicate -> sort/index -> conversion QC against lambda + pUC19 spike-ins.
# The genome FASTA and the aligner backend are the versions that matter; the index
# backend (--bowtie2 here) must match the bismark alignment backend.
#
# For RRBS: add --rrbs to trim_galore AND remove the deduplicate_bismark step
# (MspI fixed fragment ends are not PCR duplicates).
# For PBAT/scBS: add --pbat to bismark; for non-directional: add --non_directional.

set -euo pipefail

GENOME_DIR=genome
READS_DIR=fastq
OUTPUT_DIR=aligned
LAMBDA_DIR=lambda      # unmethylated lambda phage genome: residual %meth = under-conversion
PUC19_DIR=puc19        # CpG-methylated pUC19 genome: %CpG called unmethylated = over-conversion
THREADS=4              # bismark --parallel instances PER direction; total CPU scales up several-fold

mkdir -p "$OUTPUT_DIR"

bismark_genome_preparation --bowtie2 "$GENOME_DIR"

trim_galore --paired --output_dir "$READS_DIR" \
    "${READS_DIR}/sample_R1.fastq.gz" "${READS_DIR}/sample_R2.fastq.gz"

bismark --genome "$GENOME_DIR" \
    -1 "${READS_DIR}/sample_R1_val_1.fq.gz" \
    -2 "${READS_DIR}/sample_R2_val_2.fq.gz" \
    --bowtie2 \
    --parallel "$THREADS" \
    -o "$OUTPUT_DIR"

deduplicate_bismark --paired --bam \
    "${OUTPUT_DIR}/sample_R1_val_1_bismark_bt2_pe.bam"

samtools sort "${OUTPUT_DIR}/sample_R1_val_1_bismark_bt2_pe.deduplicated.bam" \
    -o "${OUTPUT_DIR}/sample.sorted.bam"
samtools index "${OUTPUT_DIR}/sample.sorted.bam"

# Conversion QC, both directions (align the same reads to each spike-in genome).
# Target: lambda non-conversion <=1% residual; pUC19 ~96-98% CpG methylated.
for spike in "$LAMBDA_DIR" "$PUC19_DIR"; do
    bismark_genome_preparation --bowtie2 "$spike"
    bismark --genome "$spike" \
        -1 "${READS_DIR}/sample_R1_val_1.fq.gz" \
        -2 "${READS_DIR}/sample_R2_val_2.fq.gz" \
        --bowtie2 -o "${spike}_qc"
done

# Mapping efficiency, per-context %methylation, and conversion lines live in the report.
cat "${OUTPUT_DIR}"/*_PE_report.txt
