#!/bin/bash
# SQANTI3 quality control for Iso-Seq isoforms

TRANSCRIPTS=$1     # clustered.hq.fasta from Iso-Seq3
ANNOTATION=$2      # Reference GTF
GENOME=$3          # Reference genome FASTA
OUTPUT_PREFIX=${4:-"sqanti"}

# Run SQANTI3 QC
# Classifies isoforms into categories:
# - FSM: Full Splice Match (known isoform)
# - ISM: Incomplete Splice Match (truncated)
# - NIC: Novel In Catalog (novel combination of known junctions)
# - NNC: Novel Not in Catalog (contains novel junctions)
# - Genic: overlaps gene but no splice match
# - Antisense, Intergenic, etc.

sqanti3_qc.py \
    $TRANSCRIPTS \
    $ANNOTATION \
    $GENOME \
    -o $OUTPUT_PREFIX \
    --report both \
    --saturation

# Key outputs:
# - ${OUTPUT_PREFIX}_classification.txt: isoform categories
# - ${OUTPUT_PREFIX}_junctions.txt: splice junction info
# - ${OUTPUT_PREFIX}_report.html: QC summary

# Filter for high-quality novel isoforms
# NIC/NNC with multiple supporting reads
awk -F'\t' '$6 == "NIC" || $6 == "NNC" {
    if ($9 >= 2) print  # min 2 full-length reads
}' ${OUTPUT_PREFIX}_classification.txt > ${OUTPUT_PREFIX}_novel_isoforms.txt

echo "Classification complete."
echo "Total isoforms: $(wc -l < ${OUTPUT_PREFIX}_classification.txt)"
echo "Novel isoforms (NIC/NNC, FL>=2): $(wc -l < ${OUTPUT_PREFIX}_novel_isoforms.txt)"

# Summary by category
echo "Category breakdown:"
cut -f6 ${OUTPUT_PREFIX}_classification.txt | sort | uniq -c | sort -rn
