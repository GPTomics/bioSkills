#!/bin/bash
# Create promoter regions from gene annotations

# Create sample data
cat > genes.bed << 'EOF'
chr1	1000	5000	GENE1	0	+
chr1	10000	15000	GENE2	0	-
chr2	5000	8000	GENE3	0	+
EOF

cat > genome.txt << 'EOF'
chr1	50000
chr2	50000
EOF

echo "=== Input Genes ==="
cat genes.bed

# Step 1: Extract TSS (transcript start sites)
echo -e "\n=== TSS Positions ==="
awk -v OFS='\t' '{
    if ($6 == "+") print $1, $2, $2+1, $4, $5, $6;
    else print $1, $3-1, $3, $4, $5, $6;
}' genes.bed > tss.bed
cat tss.bed

# Step 2: Create promoter regions (2kb upstream, 500bp downstream)
echo -e "\n=== Promoters (2kb upstream, 500bp downstream) ==="
bedtools slop -i tss.bed -g genome.txt -l 2000 -r 500 -s > promoters.bed
cat promoters.bed

# Step 3: Get upstream flanks only (2kb upstream of TSS)
echo -e "\n=== Upstream Regions Only (2kb) ==="
bedtools flank -i tss.bed -g genome.txt -l 2000 -r 0 -s > upstream.bed
cat upstream.bed

# Step 4: Get gene body + extensions
echo -e "\n=== Extended Gene Bodies (+1kb each side) ==="
bedtools slop -i genes.bed -g genome.txt -b 1000 > extended_genes.bed
cat extended_genes.bed

# Cleanup
rm -f genes.bed genome.txt tss.bed promoters.bed upstream.bed extended_genes.bed

echo -e "\n=== Commands Summary ==="
echo "Extract TSS:        awk '{if (\$6==\"+\") ... else ...}' genes.bed > tss.bed"
echo "Create promoters:   bedtools slop -i tss.bed -g genome.txt -l 2000 -r 500 -s"
echo "Upstream only:      bedtools flank -i tss.bed -g genome.txt -l 2000 -r 0 -s"
echo "Extend genes:       bedtools slop -i genes.bed -g genome.txt -b 1000"
