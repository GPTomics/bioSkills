#!/usr/bin/env bash
# Reference: DADA2 1.30+, cutadapt 4.4+ | Verify API if version differs
# OBITools3 eDNA metabarcoding pipeline: paired-end FASTQ to taxonomy table
# Marker: COI with mlCOIintF/jgHCO2198 primers
# Adjust length filters and reference database for other markers

set -euo pipefail

INPUT_R1="raw_R1.fastq.gz"
INPUT_R2="raw_R2.fastq.gz"
NGSFILTER="ngsfilter.txt"
REFDB_FASTA="MIDORI2_COI_reference.fasta"
DMS="EDNA"

# --- Step 1: Primer removal with cutadapt ---
# mlCOIintF forward / jgHCO2198 reverse
# min_overlap=20: prevents spurious matches to short random sequences
cutadapt -g 'GGWACWGGWTGAACWGTWTAYCCYCC;min_overlap=20' \
         -G 'TAIACYTCIGGRTGICCRAARAAYCA;min_overlap=20' \
         --discard-untrimmed --pair-filter=any \
         -j 4 --minimum-length 100 \
         -o trimmed_R1.fastq.gz -p trimmed_R2.fastq.gz \
         "$INPUT_R1" "$INPUT_R2"

# --- Step 2: Import into OBITools3 ---
obi import --fastq-input trimmed_R1.fastq.gz "${DMS}/reads1"
obi import --fastq-input trimmed_R2.fastq.gz "${DMS}/reads2"

# --- Step 3: Paired-end alignment ---
obi alignpairedend -R "${DMS}/reads2" "${DMS}/reads1" "${DMS}/aligned"

# score >= 50: alignment quality threshold removing poorly overlapping pairs
obi grep -p 'sequence["score"] >= 50' "${DMS}/aligned" "${DMS}/good_align"

# --- Step 4: Length filtering ---
# 200-500 bp: typical COI amplicon range for mlCOIintF/jgHCO2198 (~313 bp target)
obi grep -p 'len(sequence) >= 200 and len(sequence) <= 500' \
    "${DMS}/good_align" "${DMS}/length_ok"

# --- Step 5: Demultiplex (skip if already demultiplexed) ---
if [ -f "$NGSFILTER" ]; then
    obi ngsfilter -t "$NGSFILTER" -u "${DMS}/unassigned" \
        "${DMS}/length_ok" "${DMS}/demux"
    NEXT="${DMS}/demux"
else
    NEXT="${DMS}/length_ok"
fi

# --- Step 6: Dereplicate ---
obi uniq "${NEXT}" "${DMS}/derep"

# count >= 2: removes singletons (likely sequencing errors)
obi grep -p 'sequence["count"] >= 2' "${DMS}/derep" "${DMS}/no_singletons"

# --- Step 7: Denoise ---
# ratio=0.05: sequences with <5% abundance of a more abundant 1-mismatch parent
# are merged into the parent (removes PCR/sequencing errors)
# -s merged_sample: per-sample denoising; -r 0.05: ratio threshold; -H: head sequences only
obi clean -s merged_sample -r 0.05 -H "${DMS}/no_singletons" "${DMS}/denoised"

# --- Step 8: Import reference database ---
obi import --fasta-input "$REFDB_FASTA" "${DMS}/refdb"
obi import --taxdump taxdump.tar.gz "${DMS}/taxonomy"

# --- Step 9: Taxonomy assignment ---
# ecotag uses LCA algorithm against the reference database
obi ecotag -R "${DMS}/refdb" --taxonomy "${DMS}/taxonomy" \
    "${DMS}/denoised" "${DMS}/assigned"

# --- Step 10: Final filtering ---
# best_identity >= 0.97: species-level COI barcode gap threshold
# best_identity >= 0.95 and < 0.97: genus-level assignment
obi grep -p 'sequence["best_identity"] >= 0.95' \
    "${DMS}/assigned" "${DMS}/final"

# --- Step 11: Export ---
obi export --tab-output "${DMS}/final" > species_table.tsv

echo "Pipeline complete. Species table written to species_table.tsv"
echo "Total MOTUs:"
wc -l < species_table.tsv
