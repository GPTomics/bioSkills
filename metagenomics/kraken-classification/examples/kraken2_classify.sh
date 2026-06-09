#!/bin/bash
# Reference: Kraken2 2.1.3+, Bracken 2.9+ | Verify API if version differs
# Classify paired-end shotgun reads with the precision levers the defaults omit,
# then re-estimate species abundance with Bracken. Output is a database-conditioned
# similarity ledger: tune the database, confidence, and false-positive control.
set -euo pipefail

KRAKEN_DB="/path/to/kraken2_db"   # the database defines what can be detected
READS_R1="sample_R1.fastq.gz"
READS_R2="sample_R2.fastq.gz"
OUTPUT_DIR="kraken_output"
SAMPLE="sample"
THREADS=8
CONFIDENCE=0.1                    # raise from default 0 to suppress single-k-mer false positives (Liu 2024)
HIT_GROUPS=2                      # require >=2 distinct hit regions (raise to 3 for clinical samples)
READLEN=150                       # Bracken kmer_distrib and -r MUST match the actual read length and the DB build

mkdir -p "$OUTPUT_DIR"

kraken2 --db "$KRAKEN_DB" \
    --threads "$THREADS" \
    --paired \
    --gzip-compressed \
    --confidence "$CONFIDENCE" \
    --minimum-hit-groups "$HIT_GROUPS" \
    --use-names \
    --report "${OUTPUT_DIR}/${SAMPLE}.kreport" \
    --output "${OUTPUT_DIR}/${SAMPLE}.kraken" \
    "$READS_R1" "$READS_R2"

echo "Top 10 species by read count (NOT abundance - run Bracken, then mind genome-size bias):"
awk '$4 == "S"' "${OUTPUT_DIR}/${SAMPLE}.kreport" | sort -k1 -nr | head -10

# Re-estimate species abundance. Bracken fixes the wrong-rank problem only; it never
# removes false positives, so the confidence and hit-group filters above must run first.
bracken -d "$KRAKEN_DB" \
    -i "${OUTPUT_DIR}/${SAMPLE}.kreport" \
    -o "${OUTPUT_DIR}/${SAMPLE}.bracken" \
    -w "${OUTPUT_DIR}/${SAMPLE}.bracken.kreport" \
    -r "$READLEN" \
    -l S \
    -t 10                        # redistribution floor: drops taxa with fewer than 10 clade-level reads, strict < (not a confidence)

classified=$(grep -c "^C" "${OUTPUT_DIR}/${SAMPLE}.kraken" || true)
unclassified=$(grep -c "^U" "${OUTPUT_DIR}/${SAMPLE}.kraken" || true)
total=$((classified + unclassified))
pct=$(echo "scale=2; $classified * 100 / $total" | bc)
# A low rate can mean novel taxa OR host contamination OR the wrong database - diagnose which.
echo "Classification rate: ${pct}% (${classified}/${total})"
