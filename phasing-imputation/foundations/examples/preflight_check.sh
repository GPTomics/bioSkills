#!/bin/bash
# Reference: bcftools 1.19+ | Verify API if version differs
#
# Pipeline-readiness gate for phasing/imputation. Every stage downstream
# (phase -> impute -> filter) fails SILENTLY when the input is misaligned to the
# panel, so this script reports the conditions that must hold BEFORE phasing:
# genome-build hint, biallelic-SNP fraction, normalization, the strand-ambiguous
# (palindromic) SNP burden, and the MAF spectrum. It only READS the VCF and prints
# a report; it does not modify data, download panels, or run an engine, so it is
# safe to run on a small VCF offline.
#
# Usage: ./preflight_check.sh study.vcf.gz
set -euo pipefail

VCF=${1:?'usage: preflight_check.sh study.vcf.gz'}

PALINDROME_WARN_FRAC=0.05   # >5% A/T+C/G SNPs is a real strand-flip risk near MAF 0.5; these flip silently and the allele-frequency check cannot resolve them
RARE_MAF=0.01               # variants below this impute worst and dominate the low-INFO tail; reported so the panel/strategy choice is informed

echo "== preflight: ${VCF} =="

n_total=$(bcftools view -H "${VCF}" | wc -l | tr -d ' ')
n_samples=$(bcftools query -l "${VCF}" | wc -l | tr -d ' ')
echo "samples: ${n_samples}    records: ${n_total}"

first_chrom=$(bcftools view -H "${VCF}" | head -1 | cut -f1)
if [[ "${first_chrom}" == chr* ]]; then
    echo "chrom naming: chr-prefixed (${first_chrom}) -- panel must use the SAME convention and genome build"
else
    echo "chrom naming: numeric (${first_chrom}) -- panel must use the SAME convention and genome build"
fi

n_biallelic_snp=$(bcftools view -H -m2 -M2 -v snps "${VCF}" | wc -l | tr -d ' ')
echo "biallelic SNPs: ${n_biallelic_snp} / ${n_total} -- phasing/imputation engines want biallelic; split multiallelics with 'bcftools norm -m -any' first"

n_palindrome=$(bcftools view -H -m2 -M2 -v snps "${VCF}" \
    | awk '($4=="A"&&$5=="T")||($4=="T"&&$5=="A")||($4=="C"&&$5=="G")||($4=="G"&&$5=="C")' \
    | wc -l | tr -d ' ')
frac_palindrome=$(awk -v a="${n_palindrome}" -v b="${n_biallelic_snp}" 'BEGIN{print (b>0)? a/b : 0}')
echo "strand-ambiguous (A/T, C/G) SNPs: ${n_palindrome} (frac ${frac_palindrome})"
awk -v f="${frac_palindrome}" -v w="${PALINDROME_WARN_FRAC}" 'BEGIN{if(f>w) print "  WARN: high palindromic burden -- resolve strand by allele frequency vs the panel; drop those near MAF 0.5"}'

bcftools +fill-tags "${VCF}" -Ou -- -t MAF 2>/dev/null \
    | bcftools query -f '%INFO/MAF\n' 2>/dev/null \
    | awk -v r="${RARE_MAF}" '$1!="."{n++; if($1<r) rare++} END{if(n>0) printf "MAF spectrum: %d variants, %d (%.1f%%) below MAF %s -- the low-INFO tail; panel ancestry/size sets how well these impute\n", n, rare, 100*rare/n, r}'

echo "== gate: confirm build+strand match the panel, split multiallelics, THEN phase =="
