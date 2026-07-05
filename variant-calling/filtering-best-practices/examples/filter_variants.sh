#!/bin/bash
# Reference: GATK 4.6+, bcftools 1.19+ | Verify API if version differs
# Hard-filters a germline VCF: splits by type, applies GATK-recommended thresholds
# (SNPs and indels differ), guards RankSum annotations so hom-alt sites survive, merges.
set -euo pipefail

INPUT_VCF=$1
OUTPUT_PREFIX=${2:-filtered}

echo "=== Splitting SNPs and Indels (their error processes and thresholds differ) ==="
# This split assumes a caller (e.g. GATK HaplotypeCaller) that emits no MNP or mixed records:
# `-v snps` + `-v indels` silently DROPS MNPs and mixed SNP+indel sites. If the caller emits
# them (FreeBayes, Octopus), atomize first (bcftools norm --atomize) or add `-v mnps,other`.
bcftools view -v snps "$INPUT_VCF" -Oz -o "${OUTPUT_PREFIX}_snps_raw.vcf.gz"
bcftools view -v indels "$INPUT_VCF" -Oz -o "${OUTPUT_PREFIX}_indels_raw.vcf.gz"

# GATK-recommended lenient SNP thresholds. The "|| INFO/X = \".\"" guard on every RankSum
# term replicates GATK's "missing => PASS" rule: MQRankSum/ReadPosRankSum are undefined at
# hom-alt sites, so without the guard every hom-alt SNP would silently fail and vanish.
echo "=== Filtering SNPs ==="
bcftools filter -i '
    QUAL >= 30 &&
    (INFO/QD >= 2.0 || INFO/QD = ".") &&
    (INFO/FS <= 60.0 || INFO/FS = ".") &&
    (INFO/MQ >= 40.0 || INFO/MQ = ".") &&
    (INFO/MQRankSum >= -12.5 || INFO/MQRankSum = ".") &&
    (INFO/ReadPosRankSum >= -8.0 || INFO/ReadPosRankSum = ".") &&
    (INFO/SOR <= 3.0 || INFO/SOR = ".")
' "${OUTPUT_PREFIX}_snps_raw.vcf.gz" -Oz -o "${OUTPUT_PREFIX}_snps_filtered.vcf.gz"

# Indel thresholds differ deliberately: FS loosened 60->200 (real indels have messier local
# alignments => higher legit strand bias), ReadPosRankSum tightened -8->-20 (spurious indels
# cluster at read ends), and MQ/MQRankSum dropped (less diagnostic for indels).
echo "=== Filtering Indels ==="
bcftools filter -i '
    QUAL >= 30 &&
    (INFO/QD >= 2.0 || INFO/QD = ".") &&
    (INFO/FS <= 200.0 || INFO/FS = ".") &&
    (INFO/ReadPosRankSum >= -20.0 || INFO/ReadPosRankSum = ".") &&
    (INFO/SOR <= 10.0 || INFO/SOR = ".")
' "${OUTPUT_PREFIX}_indels_raw.vcf.gz" -Oz -o "${OUTPUT_PREFIX}_indels_filtered.vcf.gz"

echo "=== Merging filtered variants ==="
bcftools concat "${OUTPUT_PREFIX}_snps_filtered.vcf.gz" "${OUTPUT_PREFIX}_indels_filtered.vcf.gz" | \
    bcftools sort -Oz -o "${OUTPUT_PREFIX}_all_filtered.vcf.gz"
bcftools index -t "${OUTPUT_PREFIX}_all_filtered.vcf.gz"

# Ti/Tv is the fastest gross-error smell test: below ~2.0 (WGS) or ~3.0 (WES) means random
# errors (Ti/Tv ~0.5) are diluting the callset, i.e. filtering was too loose.
echo "=== Statistics (before / after) ==="
bcftools stats "$INPUT_VCF" | grep -E "^SN|^TSTV" | head -20
echo ""
bcftools stats "${OUTPUT_PREFIX}_all_filtered.vcf.gz" | grep -E "^SN|^TSTV" | head -20

TSTV=$(bcftools stats "${OUTPUT_PREFIX}_all_filtered.vcf.gz" | grep "^TSTV" | cut -f5)
echo ""
echo "Ti/Tv ratio: $TSTV (expected ~2.0-2.1 for WGS, ~3.0-3.3 for exomes)"
echo "Output: ${OUTPUT_PREFIX}_all_filtered.vcf.gz"
