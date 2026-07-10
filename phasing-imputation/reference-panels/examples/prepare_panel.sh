#!/bin/bash
# Reference: bcftools 1.19+, PLINK 1.9+ | Verify API if version differs
#
# Prepare a study VCF for imputation against a reference panel: normalize, reconcile
# chromosome naming and genome build, run the strand/allele harmonization gate, and
# build the engine-specific panel format. The DOWNLOAD and FORMAT-BUILD steps need
# multi-GB panels (and HRC/TOPMed are server-only, NOT free wget targets) and are
# shown as commented reference, not executed -- only the harmonization steps run on a
# local study VCF + a panel sites file. Outputs go to an OUT dir the caller deletes.
set -euo pipefail

STUDY=${1:?'usage: prepare_panel.sh study.vcf.gz panel_sites.vcf.gz reference.fa'}
PANEL_SITES=${2:?'panel sites/legend VCF (site-only) for the harmonization check'}
REF_FA=${3:?'reference FASTA matching the PANEL build'}
OUT=${OUT:-panel_prep_out}

# shellcheck disable=SC2034  # Rayner HRC-check reads these thresholds from its own config; surfaced here as the documented defaults
PALINDROME_MAF=0.4   # Rayner-check default: A/T and C/G SNPs above this are dropped; strand is unresolvable and frequency is too near 0.5 to disambiguate
# shellcheck disable=SC2034
AF_DIFF=0.2          # flag variants whose study-vs-panel allele-frequency gap exceeds this; a large gap signals strand, build, or ancestry mismatch

mkdir -p "${OUT}"

# Acquire the panel (reference only -- panels are multi-GB; HRC/TOPMed are server-only):
#   1000G/HGDP+1kGP: public download from IGSR/gnomAD
#   HRC:    Michigan Imputation Server (upload study, server imputes)
#   TOPMed: TOPMed/BioData Catalyst server (never downloadable)

# 1. Normalize: split multiallelics + left-align BEFORE allele harmonization, else the same indel two ways will not match the panel
bcftools norm -m -any -f "${REF_FA}" "${STUDY}" -Oz -o "${OUT}/study.norm.vcf.gz"
bcftools index "${OUT}/study.norm.vcf.gz"

# 2. Reconcile chromosome naming to the panel convention (GRCh38/TOPMed use 'chr', GRCh37 panels do not)
study_chr=$(bcftools view -H "${OUT}/study.norm.vcf.gz" | head -1 | cut -f1)
panel_chr=$(bcftools view -H "${PANEL_SITES}" | head -1 | cut -f1)
echo "study chrom: ${study_chr}    panel chrom: ${panel_chr}"
# If they differ, build a rename map and: bcftools annotate --rename-chrs map.txt

# 3. Restrict to biallelic SNPs/indels present in the panel sites (the join the engine will need)
bcftools isec -n=2 -w1 "${OUT}/study.norm.vcf.gz" "${PANEL_SITES}" -Oz -o "${OUT}/study.in_panel.vcf.gz"
bcftools index "${OUT}/study.in_panel.vcf.gz"

# 4. Strand/allele alignment to the panel reference. The field-standard gate is Will Rayner's
#    HRC-1000G-check-bim.pl, which emits Run-plink.sh (positions, ref/alt, strand; drops
#    palindromic MAF>${PALINDROME_MAF}; flags AF gap>${AF_DIFF}) and a FreqPlot. bcftools +fixref
#    is a lighter alternative that flips strand against the reference:
bcftools +fixref "${OUT}/study.in_panel.vcf.gz" -Oz -o "${OUT}/study.fixref.vcf.gz" -- \
    -f "${REF_FA}" -m flip
bcftools +fixref "${OUT}/study.fixref.vcf.gz" -- -f "${REF_FA}" -m stats

# 5. Build the engine panel format (reference -- needs the FULL haplotype panel, not sites):
#   Minimac4: minimac4 --compress-reference ref.vcf.gz > ref.msav
#   Beagle:   java -jar bref3.jar ref.vcf.gz > ref.bref3   # needs fully phased, |-separated, per chrom
#   IMPUTE5:  imp5Converter --h ref.vcf.gz --r chr20 --o ref.chr20.imp5

echo "== harmonized study written to ${OUT}/study.fixref.vcf.gz; read the AF concordance before imputing =="
