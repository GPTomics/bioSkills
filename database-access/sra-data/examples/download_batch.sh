#!/bin/bash
# Download multiple SRA runs from a file

ACCESSION_FILE="${1:-accessions.txt}"
OUTDIR="${2:-./fastq}"
THREADS="${3:-4}"

if [[ ! -f "$ACCESSION_FILE" ]]; then
    echo "Error: File $ACCESSION_FILE not found"
    echo "Usage: $0 <accessions.txt> [output_dir] [threads]"
    exit 1
fi

mkdir -p "$OUTDIR"

count=0
total=$(grep -cv '^#\|^$' "$ACCESSION_FILE")

while read -r SRR; do
    [[ -z "$SRR" || "$SRR" == \#* ]] && continue
    ((count++))

    echo "[$count/$total] Downloading $SRR..."

    fasterq-dump "$SRR" -O "$OUTDIR" -e "$THREADS" -p --skip-technical

    if [[ $? -eq 0 ]]; then
        echo "  Success: $SRR"
    else
        echo "  Failed: $SRR"
    fi
done < "$ACCESSION_FILE"

echo "Completed $count downloads"
