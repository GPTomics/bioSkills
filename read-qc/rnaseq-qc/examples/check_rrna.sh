#!/bin/bash
# Check rRNA contamination with SortMeRNA

FASTQ=$1
RRNA_DB=${2:-/path/to/sortmerna/rRNA_databases/smr_v4.3_default_db.fasta}
THREADS=${3:-8}

if [ -z "$FASTQ" ]; then
    echo "Usage: $0 <reads.fastq.gz> [rRNA_db] [threads]"
    exit 1
fi

NAME=$(basename $FASTQ .fastq.gz)
NAME=$(basename $NAME .fq.gz)

mkdir -p sortmerna_tmp

sortmerna \
    --ref $RRNA_DB \
    --reads $FASTQ \
    --aligned ${NAME}_rRNA \
    --other ${NAME}_non_rRNA \
    --fastx \
    --threads $THREADS \
    --workdir sortmerna_tmp

total=$(zcat $FASTQ | grep -c "^@" || cat $FASTQ | grep -c "^@")

if [ -f "${NAME}_rRNA.fastq" ]; then
    rrna=$(grep -c "^@" ${NAME}_rRNA.fastq)
elif [ -f "${NAME}_rRNA.fq" ]; then
    rrna=$(grep -c "^@" ${NAME}_rRNA.fq)
else
    rrna=0
fi

pct=$(echo "scale=2; $rrna / $total * 100" | bc)

echo "=== rRNA Check: $NAME ==="
echo "Total reads: $total"
echo "rRNA reads: $rrna"
echo "rRNA percentage: ${pct}%"

if (( $(echo "$pct > 20" | bc -l) )); then
    echo "WARNING: High rRNA contamination (>20%)"
elif (( $(echo "$pct > 10" | bc -l) )); then
    echo "NOTE: Moderate rRNA contamination (10-20%)"
else
    echo "OK: rRNA within acceptable range (<10%)"
fi

rm -rf sortmerna_tmp
