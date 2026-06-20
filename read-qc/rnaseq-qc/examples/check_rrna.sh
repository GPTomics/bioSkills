#!/bin/bash
# Reference: SortMeRNA 4.3+ | Verify API if version differs
# Check rRNA contamination with SortMeRNA. High rRNA = failed depletion / poly-A selection
# (a prep-efficiency readout); filter only to recover usable depth, the real fix is re-prep.

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

# Count reads as lines/4 (a quality line can start with '@', so grep '^@' overcounts).
# SortMeRNA mirrors the input compression by default, so the aligned file may be .fastq, .fastq.gz,
# .fq, or .fq.gz -- glob for whichever exists and zcat -f handles both compressed and plain.
total=$(( $(zcat -f $FASTQ | wc -l) / 4 ))

rrna_file=$(ls ${NAME}_rRNA.fastq.gz ${NAME}_rRNA.fastq ${NAME}_rRNA.fq.gz ${NAME}_rRNA.fq 2>/dev/null | head -1)
rrna=0
[ -n "$rrna_file" ] && rrna=$(( $(zcat -f "$rrna_file" | wc -l) / 4 ))

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
