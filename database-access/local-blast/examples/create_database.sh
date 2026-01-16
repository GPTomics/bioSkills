#!/bin/bash
# Create BLAST databases from FASTA files

# Create nucleotide database
echo "Creating nucleotide database..."
makeblastdb \
    -in reference.fasta \
    -dbtype nucl \
    -out ref_nucl_db \
    -title "Reference Nucleotide Database" \
    -parse_seqids

# Create protein database
echo "Creating protein database..."
makeblastdb \
    -in proteins.fasta \
    -dbtype prot \
    -out ref_prot_db \
    -title "Reference Protein Database" \
    -parse_seqids

# Verify database was created
echo -e "\nDatabase info:"
blastdbcmd -db ref_nucl_db -info
