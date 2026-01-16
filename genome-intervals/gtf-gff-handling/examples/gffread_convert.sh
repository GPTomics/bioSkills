#!/bin/bash
# Convert between GTF and GFF3 formats using gffread

# Check if gffread is available
if ! command -v gffread &> /dev/null; then
    echo "gffread not found. Install with: conda install -c bioconda gffread"
    exit 1
fi

# Example: Create a minimal GTF for testing
cat > test.gtf << 'EOF'
chr1	HAVANA	gene	11869	14409	.	+	.	gene_id "ENSG00000223972"; gene_name "DDX11L1"; gene_type "lncRNA";
chr1	HAVANA	transcript	11869	14409	.	+	.	gene_id "ENSG00000223972"; transcript_id "ENST00000456328"; gene_name "DDX11L1";
chr1	HAVANA	exon	11869	12227	.	+	.	gene_id "ENSG00000223972"; transcript_id "ENST00000456328"; exon_number "1";
chr1	HAVANA	exon	12613	12721	.	+	.	gene_id "ENSG00000223972"; transcript_id "ENST00000456328"; exon_number "2";
EOF

echo "=== Input GTF ==="
cat test.gtf

echo -e "\n=== Convert GTF to GFF3 ==="
gffread test.gtf -o test.gff3
cat test.gff3

echo -e "\n=== Convert GFF3 back to GTF ==="
gffread test.gff3 -T -o test_converted.gtf
cat test_converted.gtf

echo -e "\n=== Validate GTF ==="
gffread -E test.gtf 2>&1

# With a genome file, you can also:
# gffread -w transcripts.fa -g genome.fa annotation.gtf  # Extract transcript sequences
# gffread -x cds.fa -g genome.fa annotation.gtf          # Extract CDS sequences
# gffread -y proteins.fa -g genome.fa annotation.gtf     # Extract protein sequences

# Cleanup
rm -f test.gtf test.gff3 test_converted.gtf

echo -e "\n=== Common gffread commands ==="
echo "GTF to GFF3:           gffread annotation.gtf -o annotation.gff3"
echo "GFF3 to GTF:           gffread annotation.gff3 -T -o annotation.gtf"
echo "Extract transcripts:   gffread -w transcripts.fa -g genome.fa annotation.gtf"
echo "Extract CDS:           gffread -x cds.fa -g genome.fa annotation.gtf"
echo "Extract proteins:      gffread -y proteins.fa -g genome.fa annotation.gtf"
echo "Validate:              gffread -E annotation.gtf"
