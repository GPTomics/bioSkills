#!/bin/bash
# Reference: NCBI BLAST+ 2.15+ | Verify API if version differs
# Run BLAST with -task chosen for the question; large hitlist + post-filter avoids the max_target_seqs trap.

set -euo pipefail

PROGRAM="${1:-blastn}"        # blastn or blastp
TASK="${2:-dc-megablast}"     # megablast | dc-megablast | blastn | blastn-short | blastp | blastp-short
DB="${3:-ref_nucl_db}"
QUERY="${4:-query.fasta}"
OUT="${5:-blast_results.tsv}"
THREADS="${6:-8}"
EVALUE="${7:-1e-10}"

# qcovs = total query coverage by all HSPs of subject; qcovhsp = best-HSP coverage.
# staxids/sscinames require v5 database with taxonomy indexed.
FMT="6 qseqid sseqid pident length qcovs qcovhsp evalue bitscore staxids sscinames stitle"

echo "=== ${PROGRAM} -task ${TASK} ==="
echo "  Query: ${QUERY}"
echo "  DB:    ${DB}"
echo "  E:     ${EVALUE}"
echo "  Threads: ${THREADS}"
echo "  hitlist: 500 (post-filter to top-N by bit-score; avoids Shah 2019 trap)"
echo

${PROGRAM} \
    -query "${QUERY}" \
    -db "${DB}" \
    -out "${OUT}" \
    -outfmt "${FMT}" \
    -task "${TASK}" \
    -evalue "${EVALUE}" \
    -num_threads "${THREADS}" \
    -max_target_seqs 500 \
    -soft_masking true

echo
echo "=== Top hit per query by bit-score ==="
# Sort by qseqid then bitscore (col 8) descending; keep first row per qseqid
sort -k1,1 -k8,8gr "${OUT}" | awk -F'\t' '!seen[$1]++' > "${OUT%.tsv}.top1.tsv"
head -5 "${OUT%.tsv}.top1.tsv"

echo
echo "=== Coverage filter (qcovs >= 0.8) ==="
awk -F'\t' '$5 >= 80' "${OUT}" | head -5
