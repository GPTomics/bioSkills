#!/bin/bash
# Reference: Hyb 1.1+, bowtie2 2.5.3+, samtools 1.19+ | Verify API if version differs
# Chimeric AGO-CLIP analysis: extract miRNA-target chimeric reads from chimeric eCLIP / CLEAR-CLIP / miR-eCLIP libraries.
# Two-stage approach: chimera identification + seed-pairing validation.

TRIMMED_R1=$1                # eCLIP-style preprocessed FASTQ
MIRNA_FA=$2                  # mature miRNA sequences FASTA (mature_human_miRNA.fa from miRBase)
MRNA_FA=$3                   # target mRNA FASTA (e.g., GENCODE protein-coding transcripts)
EXPRESSED_MIRNAS=${4:-""}    # optional: BED/TSV of miRNAs > 100 TPM in matched small-RNA-seq
OUT_PREFIX=${5:-"chimera"}
THREADS=${6:-8}

# Step 1: Build combined miRNA + mRNA index for Hyb
# Hyb expects a single FASTA with miRNAs and mRNAs concatenated
cat $MIRNA_FA $MRNA_FA > ${OUT_PREFIX}_combined.fa

# Step 2: Run Hyb in multimer (chimera) mode with bowtie2 alignment
# bowtie2 is required for short miRNA sequences (21-23 nt); BLAST too stringent
if command -v hyb > /dev/null; then
    hyb \
        in=$TRIMMED_R1 \
        db=${OUT_PREFIX}_combined.fa \
        align=bowtie2 \
        type=mim
    # Output: .blast and .hyb files with miRNA-target chimera coordinates

    # Hyb output columns: chimera_id, type, miRNA_id, miRNA_align_start, ...
    # Filter: only human mRNA targets, not other miRNAs / rRNA / contamination
    awk '$5 ~ /ENST|NM_|XM_/' ${OUT_PREFIX}.hyb > ${OUT_PREFIX}_mrna_chimeras.tsv

    echo "Total chimeras: $(wc -l < ${OUT_PREFIX}.hyb)"
    echo "Human mRNA chimeras: $(wc -l < ${OUT_PREFIX}_mrna_chimeras.tsv)"
else
    echo "Hyb pipeline not on PATH; install from github.com/gkudla/hyb"
    exit 1
fi

# Step 3: Filter by miRNA expression
# Only consider chimeras involving miRNAs expressed > 100 TPM in matched small-RNA-seq
if [ -n "$EXPRESSED_MIRNAS" ]; then
    awk 'FNR==NR {expressed[$1]=1; next} $3 in expressed' \
        $EXPRESSED_MIRNAS ${OUT_PREFIX}_mrna_chimeras.tsv \
        > ${OUT_PREFIX}_expressed_chimeras.tsv
    echo "Expressed-miRNA chimeras: $(wc -l < ${OUT_PREFIX}_expressed_chimeras.tsv)"
else
    cp ${OUT_PREFIX}_mrna_chimeras.tsv ${OUT_PREFIX}_expressed_chimeras.tsv
    echo "No miRNA expression filter applied; consider matched small-RNA-seq filter"
fi

# Step 4: Per-miRNA target count
awk '{print $3}' ${OUT_PREFIX}_expressed_chimeras.tsv | sort | uniq -c | sort -rn \
    > ${OUT_PREFIX}_mirna_target_counts.tsv
echo ""
echo "Top miRNAs by target count:"
head -10 ${OUT_PREFIX}_mirna_target_counts.tsv

# Step 5 (optional): for miR-eCLIP probe-enriched library, filter for the enriched miRNA
# Example: hsa-miR-21 enrichment
# grep "hsa-miR-21" ${OUT_PREFIX}_expressed_chimeras.tsv > ${OUT_PREFIX}_mir21_targets.tsv

# Step 6: Validate canonical seed pairing in chimera target portion
# Use a Python script or custom awk for canonical 7mer-m8 / 8mer match
# (skipped for brevity; see python scripts in clip-seq repos)

echo ""
echo "Cross-validate with TargetScan conserved predictions:"
echo "  bedtools intersect -wa -wb -a chimera_targets_3utr.bed -b targetscan_conserved.bed"
echo ""
echo "Chimeric methods recover non-canonical 3'-compensatory pairing missed by seed-only methods"
