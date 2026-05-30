#!/bin/bash
# Reference: bedtools 2.31+ | Verify API if version differs
# Build strand-aware promoters and a strand-correct nearest-gene assignment.
# The point: a promoter is a TSS window you IMPOSE (not slop on a gene body),
# and "upstream/downstream" must be signed by gene strand (-D b), not coordinate (-D ref).

set -euo pipefail

PROMOTER_UP=2000    # bp upstream of TSS; common core-promoter convention, NOT a fact -- report and tune per assay
PROMOTER_DOWN=200   # bp downstream of TSS; asymmetric on purpose (+1 nucleosome / 5'UTR sit downstream)

cat > genes.bed << 'EOF'
chr1	1000	5000	GENE1	0	+
chr1	10000	15000	GENE2	0	-
chr2	300	8000	GENE3	0	+
EOF

cat > genome.txt << 'EOF'
chr1	50000
chr2	50000
EOF

cat > peaks.bed << 'EOF'
chr1	4800	4900	peakA	100	+
chr1	16000	16100	peakB	100	+
chr2	100	200	peakC	100	+
EOF

echo "=== Step 1: TSS (start for +, end-1 for -; a -strand TSS is the BED end) ==="
awk -v OFS='\t' '{ if ($6=="+") print $1,$2,$2+1,$4,$5,$6; else print $1,$3-1,$3,$4,$5,$6 }' genes.bed > tss.bed
cat tss.bed

echo "=== Step 2: Promoters = TSS -${PROMOTER_UP} / +${PROMOTER_DOWN}, strand-aware ==="
# -s makes -l mean "upstream of the gene" on BOTH strands (adds to the END on -strand features)
bedtools slop -i tss.bed -g genome.txt -s -l "$PROMOTER_UP" -r "$PROMOTER_DOWN" > promoters.bed
cat promoters.bed

echo "=== Step 3: Flag promoters CLIPPED at a chromosome end (width != requested) ==="
# slop on a 1-bp TSS gives width = up + 1 (the TSS base) + down; a shorter width means an end clip
EXPECTED=$((PROMOTER_UP + 1 + PROMOTER_DOWN))
awk -v OFS='\t' -v expw="$EXPECTED" '{ w=$3-$2; if (w<expw) print $4, "CLIPPED width="w" expected="expw }' promoters.bed

echo "=== Step 4: Nearest gene per peak, signed by GENE strand, non-overlapping, ties->first ==="
bedtools sort -i peaks.bed > peaks.sorted.bed
bedtools sort -i genes.bed > genes.sorted.bed
# -D b: sign by the gene's strand (biology); -io: closest NON-overlapping; -t first: one row per peak
# filter the -1 no-feature sentinel before any numeric use of the distance column
bedtools closest -a peaks.sorted.bed -b genes.sorted.bed -D b -io -t first | awk '$NF != -1' > nearest.bed
cat nearest.bed

rm -f genes.bed genome.txt peaks.bed tss.bed promoters.bed peaks.sorted.bed genes.sorted.bed nearest.bed
echo "=== Done ==="
