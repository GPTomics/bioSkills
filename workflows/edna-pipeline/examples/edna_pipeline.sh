#!/usr/bin/env bash
## eDNA metabarcoding pipeline using OBITools3: FASTQ to taxonomy and diversity.
## Marker: COI (Leray primers) - adjust primers and parameters for other markers.

set -euo pipefail

# --- Configuration ---
RAW_DIR="raw_reads"
TRIMMED_DIR="trimmed"
OBI_DIR="obi_workspace"
RESULTS_DIR="results"
THREADS=8

# COI Leray primers (mlCOIintF / jgHCO2198)
# Replace with appropriate primers for other markers (12S, ITS, rbcL, 18S)
FWD_PRIMER="GGWACWGGWTGAACWGTWTAYCCYCC"
REV_PRIMER="TAIACYTCIGGRTGICCRAARAAYCA"

# Reference database for ecotag (EMBL or custom formatted)
REF_DB="embl_coi_refs"

# --minimum-length 50: discard very short fragments after primer removal
MIN_LENGTH=50
# --min-score 40: paired-end alignment quality threshold
MIN_ALIGN_SCORE=40
# --min-length 100 --max-length 500: expected COI amplicon range (313bp for Leray)
# Adjust per marker: 12S MiFish ~170bp, ITS variable (200-600bp)
MIN_MERGED_LENGTH=100
MAX_MERGED_LENGTH=500
# count >=2: remove singletons (sequencing errors); increase to 5-10 for noisy data
MIN_COUNT=2
# --min-identity 0.97: species-level for COI; 0.95 for genus-level
# Lower for understudied taxa or markers with fewer references
MIN_IDENTITY=0.97

mkdir -p "$TRIMMED_DIR" "$OBI_DIR" "$RESULTS_DIR" fastqc_output

# --- Step 1: Quality assessment ---
fastqc -t "$THREADS" -o fastqc_output/ "$RAW_DIR"/*.fastq.gz
multiqc fastqc_output/ -o "$RESULTS_DIR"/multiqc_report

# --- Step 2: Primer removal with Cutadapt ---
# --discard-untrimmed: removes reads without primers (off-target amplification)
for r1 in "$RAW_DIR"/*_R1.fastq.gz; do
    sample=$(basename "$r1" _R1.fastq.gz)
    r2="${RAW_DIR}/${sample}_R2.fastq.gz"

    cutadapt \
        -g "$FWD_PRIMER" -G "$REV_PRIMER" \
        --discard-untrimmed \
        --minimum-length "$MIN_LENGTH" \
        -j "$THREADS" \
        -o "${TRIMMED_DIR}/${sample}_R1.fastq.gz" \
        -p "${TRIMMED_DIR}/${sample}_R2.fastq.gz" \
        "$r1" "$r2" \
        > "${RESULTS_DIR}/${sample}_cutadapt.log" 2>&1
done

# QC checkpoint: reads per sample after trimming
echo "=== Reads per sample after primer removal ==="
for f in "$TRIMMED_DIR"/*_R1.fastq.gz; do
    sample=$(basename "$f" _R1.fastq.gz)
    # reads per sample >1000 expected; negative controls <100
    count=$(zcat "$f" | awk 'END{print NR/4}')
    echo "$sample: $count reads"
done

# --- Step 3: Concatenate trimmed reads and import into OBITools3 ---
# OBITools3 works on merged files; concatenate per-sample outputs from cutadapt
cat "$TRIMMED_DIR"/*_R1.fastq.gz > "$TRIMMED_DIR"/all_R1.fastq.gz
cat "$TRIMMED_DIR"/*_R2.fastq.gz > "$TRIMMED_DIR"/all_R2.fastq.gz

obi import --fastq-input "$TRIMMED_DIR"/all_R1.fastq.gz "$OBI_DIR"/reads1
obi import --fastq-input "$TRIMMED_DIR"/all_R2.fastq.gz "$OBI_DIR"/reads2

# --- Step 4: Paired-end alignment ---
obi alignpairedend -R "$OBI_DIR"/reads2 "$OBI_DIR"/reads1 "$OBI_DIR"/aligned

# Filter by alignment score (removes poorly overlapping pairs)
obi grep -p "sequence[\"score\"] >= ${MIN_ALIGN_SCORE}" "$OBI_DIR"/aligned "$OBI_DIR"/score_filtered

# Filter by merged length (marker-dependent range)
obi grep -p "len(sequence) >= ${MIN_MERGED_LENGTH} and len(sequence) <= ${MAX_MERGED_LENGTH}" \
    "$OBI_DIR"/score_filtered "$OBI_DIR"/length_filtered

# --- Step 5: Dereplication ---
obi uniq "$OBI_DIR"/length_filtered "$OBI_DIR"/dereplicated

# Remove singletons and low-count sequences
obi grep -p "sequence[\"count\"] >= ${MIN_COUNT}" "$OBI_DIR"/dereplicated "$OBI_DIR"/denoised

# --- Step 6: Denoise (remove PCR/sequencing errors) ---
# ratio 0.05: sequences <5% abundance of a 1-mismatch parent are merged
obi denoise -r 0.05 "$OBI_DIR"/denoised "$OBI_DIR"/cleaned

echo "=== Sequences after chimera removal ==="
obi stat "$OBI_DIR"/cleaned

# --- Step 7: Taxonomy assignment ---
# ecotag uses LCA algorithm against reference database
obi ecotag -R "$REF_DB" "$OBI_DIR"/cleaned "$OBI_DIR"/assigned

# Filter by assignment quality
obi grep -p "sequence[\"best_identity\"] >= ${MIN_IDENTITY}" \
    "$OBI_DIR"/assigned "$OBI_DIR"/filtered_assigned

# Export to tabular format
obi export --tab-output "$OBI_DIR"/filtered_assigned > "$RESULTS_DIR"/taxonomy_table.tsv

# QC checkpoint: assignment rate
total=$(wc -l < "$RESULTS_DIR"/taxonomy_table.tsv)
# COI/BOLD: >90% to phylum expected; lower for understudied markers
assigned=$(awk -F'\t' 'NR>1 && $NF != "NA" {count++} END{print count+0}' "$RESULTS_DIR"/taxonomy_table.tsv)
echo "Taxonomy assignment rate: $assigned / $total sequences"

# --- Step 8: Generate sample x species matrix ---
# Pivot taxonomy table to community matrix for downstream R analysis
awk -F'\t' 'NR>1 {print}' "$RESULTS_DIR"/taxonomy_table.tsv \
    > "$RESULTS_DIR"/taxonomy_for_r.tsv

echo ""
echo "=== Pipeline complete ==="
echo "Taxonomy table: $RESULTS_DIR/taxonomy_table.tsv"
echo "MultiQC report: $RESULTS_DIR/multiqc_report/"
echo ""
echo "Next steps: load taxonomy_table.tsv into R for contamination filtering,"
echo "diversity analysis (iNEXT), and community comparison (vegan)."
