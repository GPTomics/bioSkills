#!/bin/bash
# Run PSI-BLAST with multiple iterations

QUERY=$1
DB=${2:-nr}
ITERATIONS=${3:-3}

if [ -z "$QUERY" ]; then
    echo "Usage: $0 <query.fasta> [database] [iterations]"
    exit 1
fi

NAME=$(basename $QUERY .fasta)

psiblast \
    -query $QUERY \
    -db $DB \
    -out ${NAME}_psiblast.txt \
    -out_pssm ${NAME}.pssm \
    -out_ascii_pssm ${NAME}_pssm.txt \
    -num_iterations $ITERATIONS \
    -inclusion_ethresh 0.001 \
    -evalue 0.01 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -num_threads 4

echo "Results: ${NAME}_psiblast.txt"
echo "PSSM: ${NAME}.pssm"
echo "Hits found: $(wc -l < ${NAME}_psiblast.txt)"
