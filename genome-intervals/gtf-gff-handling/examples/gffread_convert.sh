#!/bin/bash
# Reference: gffread 0.12+, AGAT 1.4+ | Verify API if version differs
# Convert GTF<->GFF3, extract sequences, and sanitize a malformed file with AGAT.
set -euo pipefail

if ! command -v gffread &> /dev/null; then
    echo "gffread not found. Install with: conda install -c bioconda gffread"
    exit 1
fi

cat > test.gtf << 'EOF'
chr1	HAVANA	gene	11869	14409	.	+	.	gene_id "ENSG00000223972"; gene_name "DDX11L1"; gene_biotype "lncRNA";
chr1	HAVANA	transcript	11869	14409	.	+	.	gene_id "ENSG00000223972"; transcript_id "ENST00000456328"; gene_name "DDX11L1";
chr1	HAVANA	exon	11869	12227	.	+	.	gene_id "ENSG00000223972"; transcript_id "ENST00000456328"; exon_number "1";
chr1	HAVANA	exon	12613	12721	.	+	.	gene_id "ENSG00000223972"; transcript_id "ENST00000456328"; exon_number "2";
EOF

echo "=== GTF -> GFF3 (default gffread output is GFF3) ==="
gffread test.gtf -o test.gff3
cat test.gff3

echo "=== GFF3 -> GTF2 (-T) ==="
gffread test.gff3 -T -o test_converted.gtf
cat test_converted.gtf

echo "=== Validate ==="
gffread -E test.gtf 2>&1 || true

# With a genome FASTA, gffread splices segments per transcript and respects the stop-codon convention:
#   gffread -w transcripts.fa -g genome.fa annotation.gtf   # spliced exon (mature transcript)
#   gffread -x cds.fa         -g genome.fa annotation.gtf   # spliced CDS nucleotide
#   gffread -y proteins.fa    -g genome.fa annotation.gtf   # translated-CDS protein

# When a file is malformed (no ##gff-version, missing gene/exon/UTR lines, duplicate IDs, mixed conventions),
# sanitize ONCE with AGAT rather than coding around it - it reconstructs the tree and recomputes phase:
#   agat_convert_sp_gxf2gxf.pl -g messy.gff3 -o clean.gff3
#   agat_convert_sp_gff2gtf.pl -g clean.gff3 -o clean.gtf

rm -f test.gtf test.gff3 test_converted.gtf
echo "Done."
