#!/bin/bash
# Reference: NCBI BLAST+ 2.15+ | Verify API if version differs
# Build v5 BLAST databases with parse_seqids + hash_index + (optional) per-seq taxids.

set -euo pipefail

REF_NUCL="${1:-reference.fasta}"
REF_PROT="${2:-proteins.fasta}"

echo "=== Build nucleotide DB (v5, parse_seqids, hash_index) ==="
makeblastdb \
    -in "$REF_NUCL" \
    -dbtype nucl \
    -blastdb_version 5 \
    -parse_seqids \
    -hash_index \
    -title "Reference nucleotide $(date -u +%Y-%m-%d)" \
    -out ref_nucl_db

echo
echo "=== Build protein DB (v5, parse_seqids, hash_index) ==="
makeblastdb \
    -in "$REF_PROT" \
    -dbtype prot \
    -blastdb_version 5 \
    -parse_seqids \
    -hash_index \
    -title "Reference protein $(date -u +%Y-%m-%d)" \
    -out ref_prot_db

echo
echo "=== Optional: taxid-aware DB ==="
echo "If your FASTA accessions have known taxids, build with a -taxid_map tsv (seqid<TAB>taxid):"
echo "  makeblastdb -in seqs.fasta -dbtype nucl -blastdb_version 5 \\"
echo "              -parse_seqids -taxid_map seqid_taxid.tsv -out taxid_db"
echo "This enables blastn/blastp '-taxids 9606' and '-taxidlist file.txt' filtering."

echo
echo "=== Verify ==="
blastdbcmd -db ref_nucl_db -info
echo
blastdbcmd -db ref_prot_db -info
