#!/bin/bash
# Reference: racon 1.5+, medaka 2.0+, minimap2 2.26+, merqury 1.3+, meryl 1.4+ | Verify API if version differs
# Read-type-matched ONT polishing (Racon rounds -> one medaka pass) with held-out Merqury QV
# at every step. The QV curve is the stop signal; "changes made" is NOT. Polish with ONT,
# evaluate with INDEPENDENT Illumina k-mers to avoid the circularity trap.
set -euo pipefail

DRAFT=$1            # Flye/Canu ONT draft assembly (FASTA)
ONT_READS=$2       # ONT reads used to polish (FASTQ)
EVAL_READS=$3      # INDEPENDENT eval reads, ideally Illumina (different platform than polishing)
GENOME_SIZE=$4     # expected genome size in bp, e.g. 5000000 (for Merqury best_k.sh)
MEDAKA_MODEL=$5    # MUST match basecaller+chemistry+caller+version, e.g. r1041_e82_400bps_sup_v5.0.0
OUTDIR=${6:-polish_out}

THREADS=16         # tune to host; mapping + consensus are CPU-bound
MAX_ROUNDS=4       # Racon gains decay after ~1-2; >4 just flips correct bases (Vaser 2017)
PLATEAU_DELTA=0.5  # stop when QV rises by less than 0.5 between rounds (empirical plateau)

mkdir -p "$OUTDIR"

# Build the held-out k-mer DB ONCE from independent reads. K from genome size, never hardcoded.
K=$(sh "$MERQURY"/best_k.sh "$GENOME_SIZE" | tail -n1 | awk '{print int($1+0.5)}')   # round float->int; meryl needs integer k
meryl count k="$K" output "$OUTDIR"/eval.meryl "$EVAL_READS"

measure_qv() {
    # Merqury QV of $1 against the held-out k-mer DB; echoes the QV (4th column of *.qv).
    local asm=$1 tag=$2
    ( cd "$OUTDIR" && merqury.sh eval.meryl "$asm" "qv_$tag" >/dev/null 2>&1 )
    awk '{print $4}' "$OUTDIR/qv_$tag.qv"
}

current=$DRAFT
prev_qv=$(measure_qv "$current" "input")
echo "Input QV: $prev_qv"

for i in $(seq 1 "$MAX_ROUNDS"); do
    minimap2 -t "$THREADS" -ax map-ont "$current" "$ONT_READS" > "$OUTDIR/racon_$i.sam"
    racon -t "$THREADS" "$ONT_READS" "$OUTDIR/racon_$i.sam" "$current" > "$OUTDIR/racon_$i.fasta"
    current="$OUTDIR/racon_$i.fasta"

    qv=$(measure_qv "$current" "racon$i")
    gain=$(awk -v a="$qv" -v b="$prev_qv" 'BEGIN{print a-b}')
    echo "Racon round $i: QV $qv (gain $gain)"

    # Stop at the plateau (or if QV fell - polishing made it worse).
    stop=$(awk -v g="$gain" -v d="$PLATEAU_DELTA" 'BEGIN{print (g < d) ? 1 : 0}')
    if [ "$stop" -eq 1 ]; then
        echo "QV plateaued/fell after round $i; stopping Racon"
        break
    fi
    prev_qv=$qv
done

# Exactly ONE medaka pass after Racon. The model MUST match the basecaller or the consensus
# silently degrades. Omit -m to let medaka auto-detect from the FASTQ header when possible.
medaka_consensus -i "$ONT_READS" -d "$current" -o "$OUTDIR/medaka" -t "$THREADS" -m "$MEDAKA_MODEL"
final="$OUTDIR/medaka/consensus.fasta"

final_qv=$(measure_qv "$final" "final")
echo "Final QV after medaka: $final_qv"

# Honest verdict: keep the polish only if QV actually rose vs the unpolished input.
input_qv=$(awk '{print $4}' "$OUTDIR/qv_input.qv")
better=$(awk -v f="$final_qv" -v i="$input_qv" 'BEGIN{print (f > i) ? 1 : 0}')
if [ "$better" -eq 1 ]; then
    cp "$final" "$OUTDIR/polished.fasta"
    echo "Polishing improved QV ($input_qv -> $final_qv): $OUTDIR/polished.fasta"
else
    cp "$DRAFT" "$OUTDIR/polished.fasta"
    echo "Polishing did NOT improve QV ($input_qv -> $final_qv); kept the original draft"
fi
