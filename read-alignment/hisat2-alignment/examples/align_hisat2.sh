#!/usr/bin/env bash
# Reference: hisat2 2.2+, samtools 1.19+ | Verify API if version differs
# Low-memory RNA-seq alignment with HISAT2: set strandedness; --dta only when assembling transcripts.

set -euo pipefail

R1=${1:-reads_1.fastq.gz}
R2=${2:-reads_2.fastq.gz}
INDEX=${3:-hisat2_index}         # index BASENAME, not a .ht2 file
OUTPUT=${4:-aligned.bam}
STRANDNESS=${5:-RF}              # RF = dUTP/TruSeq reverse (common); FR forward; unset for unstranded
THREADS=${6:-8}

echo "=== HISAT2 RNA-seq alignment (strandness=$STRANDNESS) ==="

# --dta is intentionally OMITTED here: it suppresses short-anchor junction reads and is for StringTie/
#   Cufflinks assembly only, not plain counting. HISAT2 gives uniques MAPQ 60 (GATK-friendly, no 255 fix).
hisat2 -p "$THREADS" -x "$INDEX" \
    --rna-strandness "$STRANDNESS" \
    --rg-id sample1 --rg SM:sample1 --rg PL:ILLUMINA \
    -1 "$R1" -2 "$R2" \
    --new-summary --summary-file "${OUTPUT%.bam}_summary.txt" | \
    samtools sort -@ "$THREADS" -o "$OUTPUT" -

samtools index "$OUTPUT"

echo "Output: $OUTPUT"
cat "${OUTPUT%.bam}_summary.txt"
