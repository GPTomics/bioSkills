#!/bin/bash
# Reference: T1K 1.0.6+, samtools 1.19+, bwa-mem 0.7.17+ | Verify CLI flags if version differs
#
# T1K HLA + KIR typing from WGS/WES (Song 2023 Genome Research).
# T1K covers class I (A, B, C), class II (DRB1, DRB3/4/5, DQA1, DQB1, DPA1, DPB1),
# AND KIR genes + KIR3DL2 Bw4/Bw6 ligand prediction in a single run.
#
# Use this as the 2024-2026 general-purpose HLA tool. Reference-bundle vintage
# is critical: update t1k-build with the current IPD-IMGT/HLA quarterly release
# (Jan/Apr/Jul/Oct) for accurate non-European typing.

set -euo pipefail

INPUT_BAM="${1:-input.bam}"
SAMPLE_NAME="${2:-sample}"
OUTPUT_DIR="${3:-t1k_output}"
THREADS="${THREADS:-8}"
HLA_INDEX="${HLA_INDEX:-hla_idx}"

HLA_REGION="chr6:28000000-34000000"
HLA_ALT_CONTIGS=(
    chr6_GL000250v2_alt chr6_GL000251v2_alt chr6_GL000252v2_alt chr6_GL000253v2_alt
    chr6_GL000254v2_alt chr6_GL000255v2_alt chr6_GL000256v2_alt
)

echo "=== T1K HLA + KIR typing ==="
mkdir -p "$OUTPUT_DIR"

if [ ! -d "$HLA_INDEX" ]; then
    echo "Building T1K HLA index against current IPD-IMGT/HLA release..."
    t1k-build.pl -o "$HLA_INDEX" --download IPD-IMGT
fi

echo "Extracting MHC reads (region + alt contigs; alt-aware critical for GRCh38)..."
samtools view -bh -@ "$THREADS" "$INPUT_BAM" "$HLA_REGION" "${HLA_ALT_CONTIGS[@]}" \
    > "$OUTPUT_DIR/hla_region.bam"

samtools sort -n -@ "$THREADS" "$OUTPUT_DIR/hla_region.bam" -o "$OUTPUT_DIR/hla_sorted.bam"
samtools fastq -@ "$THREADS" \
    -1 "$OUTPUT_DIR/hla_R1.fq" \
    -2 "$OUTPUT_DIR/hla_R2.fq" \
    -s "$OUTPUT_DIR/hla_singletons.fq" \
    -0 /dev/null \
    "$OUTPUT_DIR/hla_sorted.bam"

echo "Running T1K (class I + II + KIR)..."
t1k --preset hla \
    -1 "$OUTPUT_DIR/hla_R1.fq" \
    -2 "$OUTPUT_DIR/hla_R2.fq" \
    -f "$HLA_INDEX/hla_dna_seq.fa" \
    -o "$OUTPUT_DIR/$SAMPLE_NAME" \
    --threads "$THREADS"

GENO="$OUTPUT_DIR/${SAMPLE_NAME}_genotype.tsv"
if [ ! -f "$GENO" ]; then
    echo "Error: T1K output missing: $GENO"
    exit 1
fi

echo
echo "=== HLA + KIR Genotype ==="
column -t "$GENO"

echo
echo "=== DRB1 + DRB3/4/5 Linkage QC ==="
# DR haplotype linkage rule:
#   DR1/DR8/DR10 -> none
#   DR3/DR11/DR12/DR13/DR14 -> DRB3
#   DR4/DR7/DR9 -> DRB4
#   DR15/DR16 -> DRB5
DRB1_CALLS=$(awk -F'\t' '$1=="HLA-DRB1" {print $2; print $3}' "$GENO")
DRB3_PRESENT=$(awk -F'\t' '$1=="HLA-DRB3" && $2!="-" && $2!=""' "$GENO" | wc -l)
DRB4_PRESENT=$(awk -F'\t' '$1=="HLA-DRB4" && $2!="-" && $2!=""' "$GENO" | wc -l)
DRB5_PRESENT=$(awk -F'\t' '$1=="HLA-DRB5" && $2!="-" && $2!=""' "$GENO" | wc -l)

for allele in $DRB1_CALLS; do
    fam=$(echo "$allele" | sed 's/^HLA-DRB1\*//' | cut -d: -f1)
    case "$fam" in
        01|08|10)
            expected="none"
            ;;
        03|11|12|13|14)
            expected="DRB3"
            ;;
        04|07|09)
            expected="DRB4"
            ;;
        15|16)
            expected="DRB5"
            ;;
        *)
            expected="?"
            ;;
    esac
    echo "DRB1*$fam expects $expected"
done

echo "Observed: DRB3=$DRB3_PRESENT DRB4=$DRB4_PRESENT DRB5=$DRB5_PRESENT"

echo
echo "=== Pharmacogenomic Screen (requires 4-field specificity) ==="

check_allele() {
    local pattern=$1 drug=$2 reaction=$3
    if grep -E "\b$pattern\b" "$GENO" > /dev/null 2>&1; then
        echo "WARNING: $pattern detected -- $drug $reaction risk"
    fi
}

# Note: 4-field specificity required. B*57 family includes *57:01 (risk) and *57:03 (no risk).
check_allele "HLA-B\*57:01" "Abacavir" "hypersensitivity (Mallal 2008 NEJM)"
check_allele "HLA-B\*15:02" "Carbamazepine/oxcarbazepine" "SJS/TEN (Chung 2004 Nature; Han Chinese)"
check_allele "HLA-A\*31:01" "Carbamazepine" "DRESS/SJS (McCormack 2011 NEJM; European)"
check_allele "HLA-B\*58:01" "Allopurinol" "SJS/TEN (Hung 2005 PNAS; Han Chinese)"
check_allele "HLA-B\*13:01" "Dapsone" "DDS (Zhang 2013 NEJM; Han Chinese)"
check_allele "HLA-B\*35:02" "Minocycline" "DILI (Urban 2017; note: *35:02 NOT *35:01)"
check_allele "HLA-B\*35:01" "TMP-SMX" "DILI"
check_allele "HLA-B\*15:13" "Phenytoin" "SJS (Malaysian)"

echo
echo "T1K typing complete. Output: $GENO"
