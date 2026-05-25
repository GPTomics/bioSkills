#!/bin/bash
# Reference: NCBI BLAST+ 2.15+ | Verify API if version differs
# Quick reciprocal-best-hit ortholog candidates. For principled orthology with paralog handling, use OrthoFinder (ortholog-inference skill).

set -euo pipefail

A="${1:-species_A.fasta}"
B="${2:-species_B.fasta}"
THREADS="${3:-8}"
EVALUE="${4:-1e-5}"

echo "=== Build v5 databases with parse_seqids ==="
makeblastdb -in "${A}" -dbtype prot -blastdb_version 5 -parse_seqids -hash_index -out species_A_db
makeblastdb -in "${B}" -dbtype prot -blastdb_version 5 -parse_seqids -hash_index -out species_B_db

# hitlist_size=5 instead of 1: best-hit is determined after the search by ranking on bit-score.
# Setting max_target_seqs=1 risks the Shah 2019 early-termination bias.
FMT="6 qseqid sseqid pident length qcovs evalue bitscore"

echo
echo "=== Forward: A vs B ==="
blastp -query "${A}" -db species_B_db \
    -outfmt "${FMT}" \
    -evalue "${EVALUE}" -max_target_seqs 5 -num_threads "${THREADS}" \
    -out A_vs_B.tsv

echo
echo "=== Reverse: B vs A ==="
blastp -query "${B}" -db species_A_db \
    -outfmt "${FMT}" \
    -evalue "${EVALUE}" -max_target_seqs 5 -num_threads "${THREADS}" \
    -out B_vs_A.tsv

echo
echo "=== Best by bit-score per query, both directions ==="
# Bit-score is col 7
sort -k1,1 -k7,7gr A_vs_B.tsv | awk -F'\t' '!seen[$1]++ {print $1"\t"$2}' > A_best.tsv
sort -k1,1 -k7,7gr B_vs_A.tsv | awk -F'\t' '!seen[$1]++ {print $1"\t"$2}' > B_best.tsv

echo
echo "=== Reciprocal best hits ==="
awk -F'\t' 'NR==FNR {best_B[$1]=$2; next}
            $2 in best_B && best_B[$2]==$1 {print $1"\t"$2}' \
    A_best.tsv B_best.tsv > rbh.tsv

count=$(wc -l < rbh.tsv)
echo "Found ${count} reciprocal best hit pairs"
echo
echo "Caveat: RBH on a single forward+reverse search mis-pairs paralogs from gene duplications."
echo "For principled orthology see the ortholog-inference skill (OrthoFinder, OMA, Ensembl Compara)."
echo
echo "First 10 pairs:"
head -10 rbh.tsv | column -t
