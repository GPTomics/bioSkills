#!/usr/bin/env bash
# Reference: STAR 2.7.11+, samtools 1.19+ | Verify API if version differs
# STAR RNA-seq alignment: two-pass + GeneCounts, with the MAPQ-255 fix and strandedness detection.

set -euo pipefail

R1=${1:-reads_1.fastq.gz}
R2=${2:-reads_2.fastq.gz}
GENOME_DIR=${3:-star_index}
OUTPUT_PREFIX=${4:-star_out/sample}
THREADS=${5:-8}

mkdir -p "$(dirname "$OUTPUT_PREFIX")"

echo "=== STAR RNA-seq alignment ==="

# --outSAMmapqUnique 60: STAR gives uniques MAPQ 255 ("unavailable" in SAM); GATK drops those reads.
#   Set 60 so a downstream RNA-variant step keeps them. Harmless for plain counting.
# STAR coordinate-sorts natively -- no samtools sort needed.
# --outSAMattrRGline injects read groups (SPACE-separated tags, not bwa's -R '@RG\t...'); GATK requires them.
STAR --runThreadN "$THREADS" \
    --genomeDir "$GENOME_DIR" \
    --readFilesIn "$R1" "$R2" \
    --readFilesCommand zcat \
    --outFileNamePrefix "${OUTPUT_PREFIX}_" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --quantMode GeneCounts \
    --twopassMode Basic \
    --outSAMattrRGline ID:sample1 SM:sample1 PL:ILLUMINA LB:lib1 \
    --outSAMmapqUnique 60

samtools index "${OUTPUT_PREFIX}_Aligned.sortedByCoord.out.bam"

# Infer library strandedness from ReadsPerGene.out.tab (skip the 4 N_* summary rows; compare cols 3 vs 4).
echo "=== strandedness (use the indicated column for counting) ==="
awk 'NR>4 {f+=$3; r+=$4} END {
    printf "forward(col3)=%d  reverse(col4)=%d -> %s\n", f, r,
        (f>2*r ? "FORWARD: use col 3 / -s 1" : r>2*f ? "REVERSE: use col 4 / -s 2 (dUTP/TruSeq)" : "UNSTRANDED: use col 2 / -s 0")
}' "${OUTPUT_PREFIX}_ReadsPerGene.out.tab"

echo "BAM: ${OUTPUT_PREFIX}_Aligned.sortedByCoord.out.bam"
echo "Gene counts: ${OUTPUT_PREFIX}_ReadsPerGene.out.tab"
