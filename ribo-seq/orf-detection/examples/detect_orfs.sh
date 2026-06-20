#!/bin/bash
# Reference: RiboCode 1.2+ | Verify API if version differs
# De novo ORF detection from Ribo-seq with RiboCode (periodicity-based).

set -euo pipefail

BAM=$1           # transcriptome-aligned, P-site-calibratable Ribo-seq BAM
GTF=$2
GENOME=$3
OUTDIR=${4:-ribocode_output}
ALT_STARTS=${5:-CTG,GTG}   # near-cognate starts for uORF discovery; "" to disable

mkdir -p "$OUTDIR"

# Step 1: prepare transcript annotation
prepare_transcripts -g "$GTF" -f "$GENOME" -o "$OUTDIR/annot"

# Step 2: metaplots selects periodic read lengths AND their P-site offsets into a config.
# Read lengths come from THIS step, NOT a -l flag.
metaplots -a "$OUTDIR/annot" -r "$BAM" -o "$OUTDIR/metaplots"

# Step 3: call ORFs. -l is the longest-ORF toggle (yes/no), not read lengths.
# -A adds near-cognate starts so uORFs at CUG/GUG are found.
ALT_ARG=""
if [ -n "$ALT_STARTS" ]; then
    ALT_ARG="-A $ALT_STARTS"
fi
RiboCode -a "$OUTDIR/annot" -c "$OUTDIR/metaplots_pre_config.txt" \
    $ALT_ARG -l no -p 0.05 -o "$OUTDIR/ribocode_result"

# Summarize by ORF type. RiboCode emits:
#   annotated / uORF / dORF / Overlap_uORF / Overlap_dORF / Internal / novel
# The result table is <output_name>.txt (here ribocode_result.txt).
RESULT="$OUTDIR/ribocode_result.txt"
if [ -f "$RESULT" ]; then
    echo "ORF counts by type:"
    awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="ORF_type") c=i; next}
                {types[$c]++} END{for(t in types) print "  "t": "types[t]}' "$RESULT"
    echo "Total ORFs: $(($(wc -l < "$RESULT") - 1))"
else
    echo "No ORF_result.txt found; check the RiboCode logs."
fi
