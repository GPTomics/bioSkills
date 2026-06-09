#!/bin/bash
# Reference: AMRFinderPlus 3.12+, RGI 6+ | Verify API if version differs
# Community resistome: quantify from reads (RGI bwt) AND call presence on assembled MAGs
# (AMRFinderPlus). Output is "ARG present at abundance X", never a resistance phenotype.
set -euo pipefail

READS_R1=${1:-reads_R1.fastq.gz}
READS_R2=${2:-reads_R2.fastq.gz}
MAG=${3:-mag.fasta}
OUTDIR=${4:-amr_results}
mkdir -p "$OUTDIR"

# --- Read-based quantification (no host, no context) ---
# RGI bwt maps to CARD homolog models; it CANNOT screen point-mutation SNPs.
echo "=== RGI bwt (read-based resistome) ==="
rgi bwt -1 "$READS_R1" -2 "$READS_R2" \
    -a kma -n 8 \
    -o "${OUTDIR}/resistome" --local
# Filter the gene-level table on breadth-of-coverage before trusting any call (partial hits != genes).

# --- Contig/MAG presence calling (curated per-gene cutoffs) ---
# Run --organism ONLY on a single resolved species; omit it for mixed contigs.
echo "=== AMRFinderPlus (MAG presence) ==="
amrfinder -n "$MAG" --plus --threads 8 -o "${OUTDIR}/mag_amr.tsv"

# Summarize by drug class using the HEADER name (column order varies across versions).
echo ""
echo "ARG count by drug class (presence, not phenotype):"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Class") c=i; next} c{print $c}' \
    "${OUTDIR}/mag_amr.tsv" | sort | uniq -c | sort -rn | head -10
echo "Results in ${OUTDIR}/ - report ARG presence/abundance, not resistance."
