#!/bin/bash
# Run BLAST searches with common options

DB="${1:-ref_nucl_db}"
QUERY="${2:-query.fasta}"
OUTPUT="${3:-blast_results.tsv}"
THREADS="${4:-4}"
EVALUE="${5:-1e-5}"

echo "Running BLASTN..."
echo "  Query: $QUERY"
echo "  Database: $DB"
echo "  E-value: $EVALUE"
echo "  Threads: $THREADS"

blastn \
    -query "$QUERY" \
    -db "$DB" \
    -out "$OUTPUT" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
    -evalue "$EVALUE" \
    -num_threads "$THREADS" \
    -max_target_seqs 10 \
    -max_hsps 1

echo -e "\nResults saved to: $OUTPUT"
echo "Top hits:"
head -5 "$OUTPUT" | column -t
