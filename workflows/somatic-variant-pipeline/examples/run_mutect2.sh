#!/bin/bash
# Reference: GATK 4.5+, bcftools 1.19+ | Verify API if version differs
# End-to-end Mutect2 tumor-normal somatic pipeline: call -> orientation model
# -> contamination -> FilterMutectCalls -> PASS. Template with placeholder inputs.

set -euo pipefail

TUMOR_BAM=${1:-tumor.bam}
NORMAL_BAM=${2:-normal.bam}
REFERENCE=${3:-reference.fa}
OUTPUT_PREFIX=${4:-somatic}
GNOMAD=${5:-af-only-gnomad.vcf.gz}       # AF-only gnomAD: germline prior for Mutect2
PON=${6:-pon.vcf.gz}                       # panel of normals: recurrent-artifact removal
# GetPileupSummaries needs COMMON biallelic-SNP sites (with population AF), NOT the
# af-only-gnomad used as the germline prior; the two resources serve different roles.
COMMON_SNPS=${7:-small_exac_common_3.vcf.gz}

echo "=== Somatic Variant Calling Pipeline ==="
echo "Tumor: $TUMOR_BAM | Normal: $NORMAL_BAM"

# -normal takes the read-group SM name, not the filename. awk (not sed) for portable tab
# handling: BSD/macOS sed does not interpret \t in a bracket expression.
NORMAL_NAME=$(samtools view -H "$NORMAL_BAM" | awk '/^@RG/{for(i=1;i<=NF;i++) if($i ~ /^SM:/){sub(/^SM:/,"",$i); print $i; exit}}')
echo "Normal sample name: $NORMAL_NAME"

# Step 1: somatic call (tumor + normal in one command). --f1r2-tar-gz is required
# for the orientation-bias model in Step 2; without it Step 2 is impossible.
echo "=== Step 1: Mutect2 ==="
gatk Mutect2 \
    -R "$REFERENCE" \
    -I "$TUMOR_BAM" -I "$NORMAL_BAM" -normal "$NORMAL_NAME" \
    --germline-resource "$GNOMAD" \
    --panel-of-normals "$PON" \
    --f1r2-tar-gz "${OUTPUT_PREFIX}_f1r2.tar.gz" \
    -O "${OUTPUT_PREFIX}_raw.vcf.gz"

# Step 2: orientation-bias model (oxoG C>A/G>T, FFPE C>T/G>A strand artifacts)
echo "=== Step 2: LearnReadOrientationModel ==="
gatk LearnReadOrientationModel \
    -I "${OUTPUT_PREFIX}_f1r2.tar.gz" \
    -O "${OUTPUT_PREFIX}_read_orientation.tar.gz"

# Step 3: contamination on common biallelic SNPs
echo "=== Step 3: contamination ==="
gatk GetPileupSummaries -I "$TUMOR_BAM" \
    -V "$COMMON_SNPS" -L "$COMMON_SNPS" -O "${OUTPUT_PREFIX}_tumor_pileups.table"
gatk GetPileupSummaries -I "$NORMAL_BAM" \
    -V "$COMMON_SNPS" -L "$COMMON_SNPS" -O "${OUTPUT_PREFIX}_normal_pileups.table"
gatk CalculateContamination \
    -I "${OUTPUT_PREFIX}_tumor_pileups.table" \
    -matched "${OUTPUT_PREFIX}_normal_pileups.table" \
    -O "${OUTPUT_PREFIX}_contamination.table" \
    --tumor-segmentation "${OUTPUT_PREFIX}_segments.table"

# Step 4: joint filter (consumes contamination + segmentation + orientation priors)
echo "=== Step 4: FilterMutectCalls ==="
gatk FilterMutectCalls \
    -R "$REFERENCE" \
    -V "${OUTPUT_PREFIX}_raw.vcf.gz" \
    --contamination-table "${OUTPUT_PREFIX}_contamination.table" \
    --tumor-segmentation "${OUTPUT_PREFIX}_segments.table" \
    --ob-priors "${OUTPUT_PREFIX}_read_orientation.tar.gz" \
    -O "${OUTPUT_PREFIX}_filtered.vcf.gz"

# Extract PASS somatic variants
bcftools view -f PASS "${OUTPUT_PREFIX}_filtered.vcf.gz" -Oz -o "${OUTPUT_PREFIX}_pass.vcf.gz"
bcftools index -t "${OUTPUT_PREFIX}_pass.vcf.gz"

echo "=== Summary ==="
bcftools stats "${OUTPUT_PREFIX}_pass.vcf.gz" | grep -E "^SN"
echo "PASS somatic variants: ${OUTPUT_PREFIX}_pass.vcf.gz"
