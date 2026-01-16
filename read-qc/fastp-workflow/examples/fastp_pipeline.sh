#!/bin/bash
# Complete fastp preprocessing pipeline
# Usage: ./fastp_pipeline.sh <R1.fastq.gz> <R2.fastq.gz> <output_prefix>

R1="${1}"
R2="${2}"
PREFIX="${3:-processed}"
THREADS="${4:-8}"

if [[ -z "$R1" ]] || [[ -z "$R2" ]]; then
    echo "Usage: $0 <R1.fastq.gz> <R2.fastq.gz> [output_prefix] [threads]"
    exit 1
fi

echo "=== fastp Preprocessing Pipeline ==="
echo "Input R1: $R1"
echo "Input R2: $R2"
echo "Output prefix: $PREFIX"
echo "Threads: $THREADS"

fastp \
    -i "$R1" \
    -I "$R2" \
    -o "${PREFIX}_R1.fastq.gz" \
    -O "${PREFIX}_R2.fastq.gz" \
    --detect_adapter_for_pe \
    --cut_right \
    --cut_right_window_size 4 \
    --cut_right_mean_quality 20 \
    -q 20 \
    -l 36 \
    --thread "$THREADS" \
    -h "${PREFIX}_fastp.html" \
    -j "${PREFIX}_fastp.json"

if [[ $? -eq 0 ]]; then
    echo -e "\n=== Output Files ==="
    ls -lh "${PREFIX}"*

    echo -e "\n=== Summary (from JSON) ==="
    python3 -c "
import json
with open('${PREFIX}_fastp.json') as f:
    r = json.load(f)
s = r['summary']
b = s['before_filtering']
a = s['after_filtering']
print(f\"Before: {b['total_reads']:,} reads, {b['total_bases']:,} bases\")
print(f\"After:  {a['total_reads']:,} reads, {a['total_bases']:,} bases\")
print(f\"Passed: {a['total_reads']/b['total_reads']*100:.1f}%\")
print(f\"Q20 rate: {a['q20_rate']*100:.1f}%\")
print(f\"Q30 rate: {a['q30_rate']*100:.1f}%\")
"

    echo -e "\nHTML report: ${PREFIX}_fastp.html"
else
    echo "fastp failed"
    exit 1
fi
