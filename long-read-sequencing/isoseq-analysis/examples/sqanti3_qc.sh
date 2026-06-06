#!/bin/bash
# Reference: SQANTI3 5.2+, minimap2 2.28+ | Verify API if version differs
# SQANTI3 classification + artifact filter for any long-read transcriptome (PacBio or ONT),
# with orthogonal support. The structural_category field holds LONG strings
# (full-splice_match, novel_in_catalog, ...), not short codes like FSM/NIC.
set -euo pipefail

ISOFORMS=${1:?Usage: $0 <isoforms.gff> <annotation.gtf> <genome.fa> [out_dir] [cage.bed] [polyA.txt] [sr_fofn]}
ANNOTATION=${2:?annotation GTF required}
GENOME=${3:?genome FASTA required}
OUTDIR=${4:-sqanti_out}
CAGE=${5:-}                # CAGE refTSS BED: 5' TSS validation (catches 5' degradation -> ISM)
POLYA=${6:-}               # poly-A motif list: 3' validation (catches intra-priming)
SR_FOFN=${7:-}             # short-read STAR SJ fofn: junction validation (catches RT-switch/NNC)

OPT=""
[ -n "$CAGE" ]   && OPT="$OPT --CAGE_peak $CAGE"
[ -n "$POLYA" ]  && OPT="$OPT --polyA_motif_list $POLYA"
[ -n "$SR_FOFN" ] && OPT="$OPT --short_reads $SR_FOFN"

mkdir -p "$OUTDIR"
# The isoforms positional defaults to GTF/GFF in SQANTI3; add --fasta only for a FASTA input.
sqanti3_qc.py "$ISOFORMS" "$ANNOTATION" "$GENOME" -d "$OUTDIR" -o sqanti --report both $OPT

CLASS="$OUTDIR/sqanti_classification.txt"

# Apply the artifact filter (rules = transparent thresholds; ml = random forest). The filter,
# not the discovery, is the analysis - it removes RT-switching/intra-priming/unsupported junctions.
# --gtf is the GTF/GFF input for the filter (its --isoforms flag expects FASTA/FASTQ).
sqanti3_filter.py rules "$CLASS" --gtf "$ISOFORMS" -d "$OUTDIR"

echo 'Structural-category breakdown (NIC > NNC in trust; high ISM = RNA degradation):'
col=$(head -1 "$CLASS" | tr '\t' '\n' | grep -n '^structural_category$' | cut -d: -f1)
cut -f"$col" "$CLASS" | tail -n +2 | sort | uniq -c | sort -rn
