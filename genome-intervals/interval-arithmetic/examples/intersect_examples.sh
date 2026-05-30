#!/bin/bash
# Reference: bedtools 2.31+ | Verify API if version differs
# Core interval set operations and the silent-precondition footguns.

set -euo pipefail

MERGE_DIST=100   # merge replicate peaks within 100 bp into one consensus interval (replicate-merge convention)

cat > peaks.bed << 'EOF'
chr1	100	200	peak1	100	+
chr1	300	400	peak2	200	+
chr1	500	600	peak3	150	+
chr2	100	200	peak4	250	-
EOF

cat > genes.bed << 'EOF'
chr1	150	350	geneA	0	+
chr1	550	700	geneB	0	-
chr2	50	150	geneC	0	+
EOF

cat > genome.txt << 'EOF'
chr1	248956422
chr2	242193529
EOF

echo "=== intersect -u (whole A once if it overlaps any B) ==="
bedtools intersect -a peaks.bed -b genes.bed -u

echo "=== intersect -v (A features with NO overlap) ==="
bedtools intersect -a peaks.bed -b genes.bed -v

echo "=== intersect -c (per-A count of B hits) ==="
bedtools intersect -a peaks.bed -b genes.bed -c

echo "=== intersect -wa -wb (join: whole A + whole B per pair) ==="
bedtools intersect -a peaks.bed -b genes.bed -wa -wb

echo "=== intersect -loj (left outer join: every A, NULL B if none) ==="
bedtools intersect -a peaks.bed -b genes.bed -loj

echo "=== A-vs-B fraction asymmetry: -f thresholds fraction of A ==="
bedtools intersect -a peaks.bed -b genes.bed -f 0.5 -u

echo "=== subtract (clip B out of A) vs subtract -A (drop whole A) ==="
bedtools subtract -a peaks.bed -b genes.bed
bedtools subtract -a peaks.bed -b genes.bed -A

echo "=== merge REQUIRES prior sort ==="
bedtools sort -i peaks.bed | bedtools merge -d "$MERGE_DIST" -c 4,5 -o distinct,sum

echo "=== complement: the gaps (genome file required) ==="
bedtools sort -i peaks.bed | bedtools complement -i - -g genome.txt

echo "=== -sorted is an unchecked promise: always pair with -g so order mismatch crashes ==="
bedtools intersect -a <(bedtools sort -i peaks.bed) -b <(bedtools sort -i genes.bed) -sorted -g genome.txt -u

echo "=== consensus peakset from replicates: cat | sort | merge ==="
cat peaks.bed peaks.bed | bedtools sort | bedtools merge -d "$MERGE_DIST"

echo "=== map: transfer/aggregate a B column onto each sorted A interval ==="
bedtools map -a <(bedtools sort -i genes.bed) -b <(bedtools sort -i peaks.bed) -c 5 -o mean

rm -f peaks.bed genes.bed genome.txt
