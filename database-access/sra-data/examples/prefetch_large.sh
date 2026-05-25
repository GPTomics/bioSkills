#!/bin/bash
# Reference: sra-tools 3.0+, aws-cli 2+ | Verify API if version differs
# Cloud-native SRA pull via AWS STRIDES (zero egress from EC2 in us-east-1) + 10x technical-reads support.

set -euo pipefail

SRR="${1:-SRR12345678}"
OUT="${2:-./fastq}"
THREADS="${3:-8}"
TENX="${4:-no}"  # 'yes' to keep technical reads (barcode/UMI/index) for 10x records

mkdir -p "${OUT}"

# Check if STRIDES (AWS Open Data) has the file
echo "=== Check STRIDES availability ==="
if aws s3 ls "s3://sra-pub-run-odp/sra/${SRR}/" --no-sign-request 2>/dev/null | grep -q "${SRR}"; then
    echo "  Available on AWS Open Data; pulling (free within us-east-1)"
    # STRIDES objects are unsuffixed (just SRR12345678); rename on copy.
    aws s3 cp "s3://sra-pub-run-odp/sra/${SRR}/${SRR}" "./${SRR}.sra" --no-sign-request
    SRA_PATH="./${SRR}.sra"
else
    echo "  Not on AWS; falling back to NCBI prefetch"
    prefetch "${SRR}" --max-size 200G -p
    SRA_PATH="${SRR}"
fi

echo
echo "=== Validate ==="
vdb-validate "${SRA_PATH}" || { echo "VALIDATION FAILED"; exit 1; }

echo
echo "=== fasterq-dump ==="
TECH_FLAG="--skip-technical"
if [ "${TENX}" = "yes" ]; then
    TECH_FLAG="--include-technical"
    echo "  10x mode: keeping technical reads (R1=barcode+UMI, R2=cDNA, I1=index for 10x v3)"
fi

fasterq-dump "${SRA_PATH}" \
    -O "${OUT}" \
    -e "${THREADS}" \
    -p \
    --split-files \
    ${TECH_FLAG}

echo
echo "=== Compress (fasterq-dump does NOT compress) ==="
pigz -p "${THREADS}" "${OUT}/${SRR}"_*.fastq

echo
echo "=== Cleanup ==="
rm -f "./${SRR}.sra"
rm -rf "${SRR}"  # SRA cache directory if prefetch used

echo
echo "Files:"
ls -lh "${OUT}/${SRR}"_*.fastq.gz
