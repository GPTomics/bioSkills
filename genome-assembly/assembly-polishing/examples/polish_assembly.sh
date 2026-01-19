#!/bin/bash
# Polish long-read assembly with Pilon
set -euo pipefail

ASSEMBLY=$1
R1=$2
R2=$3
OUTPUT=${4:-polished.fasta}
ROUNDS=${5:-2}

echo "=== Assembly Polishing with Pilon ==="

current=$ASSEMBLY
workdir=$(mktemp -d)

for i in $(seq 1 $ROUNDS); do
    echo "Round $i of $ROUNDS..."

    # Index and align
    bwa index $current 2>/dev/null
    bwa mem -t 16 $current $R1 $R2 2>/dev/null | samtools sort -o ${workdir}/round${i}.bam
    samtools index ${workdir}/round${i}.bam

    # Run Pilon
    pilon --genome $current \
          --frags ${workdir}/round${i}.bam \
          --output ${workdir}/pilon_round${i} \
          --changes

    current=${workdir}/pilon_round${i}.fasta

    # Report changes
    changes=$(wc -l < ${workdir}/pilon_round${i}.changes || echo "0")
    echo "  Changes made: $changes"

    if [ "$changes" -eq 0 ]; then
        echo "No more changes, stopping early"
        break
    fi
done

cp $current $OUTPUT
rm -rf $workdir

echo ""
echo "Polished assembly: $OUTPUT"
