#!/bin/bash
# Find reciprocal best BLAST hits between two species

SPECIES_A="species_A.fasta"
SPECIES_B="species_B.fasta"
THREADS="${1:-4}"
EVALUE="${2:-1e-5}"

echo "=== Creating BLAST databases ==="
makeblastdb -in "$SPECIES_A" -dbtype prot -out species_A_db -parse_seqids
makeblastdb -in "$SPECIES_B" -dbtype prot -out species_B_db -parse_seqids

echo -e "\n=== Forward BLAST: A vs B ==="
blastp -query "$SPECIES_A" -db species_B_db \
    -outfmt "6 qseqid sseqid pident evalue bitscore" \
    -evalue "$EVALUE" -max_target_seqs 1 -num_threads "$THREADS" \
    -out A_vs_B.tsv

echo "=== Reverse BLAST: B vs A ==="
blastp -query "$SPECIES_B" -db species_A_db \
    -outfmt "6 qseqid sseqid pident evalue bitscore" \
    -evalue "$EVALUE" -max_target_seqs 1 -num_threads "$THREADS" \
    -out B_vs_A.tsv

echo -e "\n=== Finding Reciprocal Best Hits ==="
# Join results where A's best hit in B has A as its best hit
awk 'NR==FNR {best_B[$1]=$2; next}
     $2 in best_B && best_B[$2]==$1 {print $1, $2, $3, $4}' \
    A_vs_B.tsv B_vs_A.tsv > reciprocal_best_hits.tsv

count=$(wc -l < reciprocal_best_hits.tsv)
echo "Found $count reciprocal best hit pairs"
echo -e "\nFirst 10 pairs:"
head -10 reciprocal_best_hits.tsv | column -t
