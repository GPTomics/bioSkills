#!/bin/bash
# Base Editing Analysis with CRISPResso2
# Demonstrates analysis of CBE and ABE experiments

# =============================================================================
# Cytosine Base Editor (CBE) Analysis
# =============================================================================

# CBE converts C->T (or G->A on opposite strand)
# Target: C at position 6 of protospacer

AMPLICON="ATGCGATCGATCGATCGATCGATCGATCGATCG"
GUIDE="TCGATCGATCGATCGAT"  # 17bp spacer
EDITED_AMPLICON="ATGCGATCGATCGTTCGATCGATCGATCGATCG"  # C->T at position 6

# Run CRISPResso2 for base editing
CRISPResso --fastq_r1 cbe_sample_R1.fastq.gz \
    --fastq_r2 cbe_sample_R2.fastq.gz \
    --amplicon_seq "$AMPLICON" \
    --guide_seq "$GUIDE" \
    --expected_hdr_amplicon_seq "$EDITED_AMPLICON" \
    --base_editor_output \
    --conversion_nuc_from C \
    --conversion_nuc_to T \
    --quantification_window_size 10 \
    --quantification_window_center -10 \
    -o cbe_results/ \
    -n cbe_sample

# Key output files:
# - Quantification_window_nucleotide_percentage_table.txt (per-position edits)
# - CRISPResso2_info.json (summary statistics)
# - Alleles_frequency_table.txt (all observed alleles)

# =============================================================================
# Adenine Base Editor (ABE) Analysis
# =============================================================================

# ABE converts A->G (or T->C on opposite strand)
ABE_AMPLICON="ATGCGATCGATCGATCGATCGATCGATCGATCG"
ABE_GUIDE="TCGATCGATCGATCGAT"
ABE_EDITED="ATGCGATCGATCGGTCGATCGATCGATCGATCG"  # A->G

CRISPResso --fastq_r1 abe_sample_R1.fastq.gz \
    --fastq_r2 abe_sample_R2.fastq.gz \
    --amplicon_seq "$ABE_AMPLICON" \
    --guide_seq "$ABE_GUIDE" \
    --expected_hdr_amplicon_seq "$ABE_EDITED" \
    --base_editor_output \
    --conversion_nuc_from A \
    --conversion_nuc_to G \
    -o abe_results/ \
    -n abe_sample

# =============================================================================
# Prime Editing Analysis
# =============================================================================

# Prime editing with pegRNA
PE_AMPLICON="ATGCGATCGATCGATCGATCGATCGATCGATCG"
PE_SPACER="TCGATCGATCGATCGAT"
PE_EXTENSION="GCTAGCTAGCTA"  # RT template + PBS
PE_EDITED="ATGCGATCGATCGATGCTAGCGATCGATCGATCG"

CRISPResso --fastq_r1 pe_sample_R1.fastq.gz \
    --fastq_r2 pe_sample_R2.fastq.gz \
    --amplicon_seq "$PE_AMPLICON" \
    --guide_seq "$PE_SPACER" \
    --expected_hdr_amplicon_seq "$PE_EDITED" \
    --prime_editing_pegRNA_extension_seq "$PE_EXTENSION" \
    -o pe_results/ \
    -n pe_sample

# =============================================================================
# Batch Processing Multiple Samples
# =============================================================================

# Process multiple samples with same amplicon/guide
for sample in sample1 sample2 sample3; do
    CRISPResso --fastq_r1 "${sample}_R1.fastq.gz" \
        --fastq_r2 "${sample}_R2.fastq.gz" \
        --amplicon_seq "$AMPLICON" \
        --guide_seq "$GUIDE" \
        --base_editor_output \
        --conversion_nuc_from C \
        --conversion_nuc_to T \
        -o "results/${sample}" \
        -n "$sample"
done

# =============================================================================
# Extract Key Metrics
# =============================================================================

echo "=== Base Editing Summary ==="
for dir in results/*/; do
    sample=$(basename "$dir")
    # Extract editing efficiency from JSON
    if [ -f "${dir}CRISPResso2_info.json" ]; then
        echo "Sample: $sample"
        python3 -c "
import json
with open('${dir}CRISPResso2_info.json') as f:
    data = json.load(f)
    print(f\"  Total reads: {data.get('aln_stats', {}).get('n_total', 'N/A')}\")
    print(f\"  Aligned: {data.get('aln_stats', {}).get('n_aligned', 'N/A')}\")
"
    fi
done

# Quality thresholds for base editing:
# - Editing efficiency >30%: Good
# - Indel rate <5%: Ideal for base editors
# - Bystander rate <10%: Often acceptable
