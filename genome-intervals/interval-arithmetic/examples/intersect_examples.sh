#!/bin/bash
# Common bedtools intersect patterns

# Setup: create sample files
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

echo "=== Input Files ==="
echo "Peaks:"
cat peaks.bed
echo -e "\nGenes:"
cat genes.bed

echo -e "\n=== Basic Intersect (overlapping portions) ==="
bedtools intersect -a peaks.bed -b genes.bed

echo -e "\n=== Peaks overlapping genes (-u) ==="
bedtools intersect -a peaks.bed -b genes.bed -u

echo -e "\n=== Peaks NOT overlapping genes (-v) ==="
bedtools intersect -a peaks.bed -b genes.bed -v

echo -e "\n=== With gene annotation (-wa -wb) ==="
bedtools intersect -a peaks.bed -b genes.bed -wa -wb

echo -e "\n=== Count overlaps (-c) ==="
bedtools intersect -a peaks.bed -b genes.bed -c

echo -e "\n=== Require 50% overlap (-f 0.5) ==="
bedtools intersect -a peaks.bed -b genes.bed -f 0.5 -u

# Cleanup
rm -f peaks.bed genes.bed
