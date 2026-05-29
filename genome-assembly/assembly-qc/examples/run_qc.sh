#!/bin/bash
# Reference: QUAST 5.2+, BUSCO 5.5+/6.x, compleasm 0.2.6+, Merqury 1.3+, meryl 1.4+, minimap2 2.26+ | Verify API if version differs
# Three-axis reference-free assembly QC: contiguity + completeness + correctness.
# No single number is quality; this runs all three axes and never computes QV on the polishing reads.
set -euo pipefail

ASSEMBLY=$1            # the assembly FASTA to evaluate
READS=$2              # ACCURATE reads (HiFi/Illumina) for the k-mer DB; ideally NOT the polishing reads
GENOME_SIZE=$3        # expected genome size in INTEGER bp (e.g. 1200000000) from GenomeScope2/flow cytometry; calN50/QUAST/best_k all want a number
LINEAGE=${4:-eukaryota}        # deepest applicable BUSCO/compleasm clade; eukaryota is a coarse fallback
OUTDIR=${5:-assembly_qc}
THREADS=${6:-16}      # tune to host; QC mapping/k-mer counting is CPU-bound

mkdir -p "$OUTDIR"

echo "=== Axis A: contiguity (auN/NGx, not bare N50) ==="
# calN50.js ships with minimap2; -L sets the genome size so NG50/auNG are genome-size-normalized
k8 "$(dirname "$(command -v minimap2)")/calN50.js" -L "$GENOME_SIZE" "$ASSEMBLY" > "$OUTDIR/calN50.txt" || \
  echo "calN50.js not found; it ships with minimap2 (k8 + calN50.js) - install minimap2, or read NG50 from QUAST below"
# --est-ref-size gives QUAST the genome size it needs for NG50/NGx (else QUAST reports N50/L50 only)
quast.py "$ASSEMBLY" --large --eukaryote --est-ref-size "$GENOME_SIZE" -t "$THREADS" -o "$OUTDIR/quast"

echo "=== Axis B: completeness (gene-space + whole-genome k-mer) ==="
# compleasm is preferred on good genomes (BUSCO under-reports via its own predictor)
compleasm run -a "$ASSEMBLY" -l "$LINEAGE" -o "$OUTDIR/compleasm" -t "$THREADS" || \
  busco -i "$ASSEMBLY" -m genome -l "${LINEAGE}_odb10" -o busco_run -c "$THREADS"

echo "=== Axis C: correctness (reference-free QV + k-mer completeness) ==="
# k from best_k.sh, NOT hardcoded; wrong k silently degrades QV/completeness
K=$(sh "$MERQURY"/best_k.sh "$GENOME_SIZE" | tail -n1 | awk '{print int($1+0.5)}')   # round float->int; meryl needs integer k
echo "best_k.sh recommends k=$K for genome size $GENOME_SIZE"
meryl count k="$K" "$READS" output "$OUTDIR/reads.meryl"
( cd "$OUTDIR" && merqury.sh reads.meryl "$(readlink -f "$ASSEMBLY")" merqury )
# merqury writes merqury.qv (overall + per-scaffold), merqury.completeness.stats, and spectra-cn plots

echo ""
echo "=== Reports in $OUTDIR ==="
echo "Contiguity : $OUTDIR/calN50.txt , $OUTDIR/quast/report.txt"
echo "Completeness: $OUTDIR/compleasm/summary.txt (or busco short_summary)"
echo "Correctness : $OUTDIR/merqury.qv , $OUTDIR/merqury.completeness.stats , $OUTDIR/merqury*.spectra-cn.*"
echo ""
echo "Reminder: a QV computed on the reads used for polishing is circular. Report all three axes."
