#!/bin/bash
# Reference: Beagle 5.4 (22Jul22), Minimac4 4.1+, bcftools 1.19+ | Verify API if version differs
#
# Array imputation against a prepared reference panel, per chromosome, emitting
# dosages with a per-variant quality. Two engines are shown: Beagle (phases the
# unphased target itself) and Minimac4 (needs a PHASED target; positional-arg CLI).
# NOT a spot-run target: it needs the engine binaries, a multi-GB panel in the
# engine format, and minute-to-hour jobs. Outputs go to an OUT dir the caller
# deletes. HRC/TOPMed panels are server-only -- see the imputation-server note.
set -euo pipefail

TARGET=${1:?'usage: run_beagle_imputation.sh study.vcf.gz panel.bref3 plink.chr20.GRCh38.map chr20'}
PANEL=${2:?'reference panel in engine format (bref3 for Beagle, msav for Minimac4)'}
MAP=${3:?'PLINK genetic map matching the panel build'}
CHR=${4:?'chromosome, e.g. chr20'}
OUT=${OUT:-imputation_out}
XMX=${XMX:-50g}      # Beagle JVM heap; scale with sample x panel size; OOM if too low. Impute per chromosome to bound memory
THREADS=${THREADS:-8}

mkdir -p "${OUT}"

# Beagle: phases the unphased target internally, then imputes. gp/ap emit GP and per-haplotype dosage
java -Xmx"${XMX}" -jar beagle.jar \
    gt="${TARGET}" \
    ref="${PANEL}" \
    map="${MAP}" \
    out="${OUT}/imputed.${CHR}" \
    nthreads="${THREADS}" \
    gp=true \
    ap=true
bcftools index "${OUT}/imputed.${CHR}.vcf.gz"

# Minimac4 alternative (the target MUST be pre-phased; positional args, NOT --refHaps):
#   minimac4 panel.msav target.phased.${CHR}.vcf.gz -o ${OUT}/imputed.${CHR}.vcf.gz -f GT,DS,HDS,GP -t ${THREADS}

# Extract the dosage-quality field (Beagle DR2; Minimac R2) for downstream filtering -> imputation-qc
bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/DR2\t%INFO/AF\n' "${OUT}/imputed.${CHR}.vcf.gz" > "${OUT}/dr2.${CHR}.txt"

# For HRC/TOPMed (server-only panels), upload a per-chromosome VCF to the Michigan or TOPMed
# Imputation Server instead of running locally; it runs Eagle2 + Minimac4 and returns GT,DS,GP + an info file.

echo "== imputed dosages in ${OUT}/imputed.${CHR}.vcf.gz; carry DS (not hard GT) to the GWAS, after filtering on DR2 =="
