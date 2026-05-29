#!/bin/bash
# Reference: SPAdes 4.0+, KMC 3.2+, GenomeScope2 2.0+ | Verify API if version differs
# Profile-first bacterial isolate assembly from Illumina paired-end reads.
# Profiling sets the size/het/repeat ceiling BEFORE assembly; SPAdes --isolate
# (not --careful, and never both) is the modern isolate default.
set -euo pipefail

R1=$1
R2=$2
OUTDIR=$3
THREADS=${4:-16}

KMER=21          # GenomeScope/KMC convention: long enough to be mostly unique, short enough for k-mer coverage
PLOIDY=2         # set 1 for a haploid bacterial isolate; 2 only if a heterozygous diploid is suspected
MEM_GB=64        # SPAdes peak RAM cap (GB); tune to host. Assembly is RAM-bound
MAX_KMER_COUNT=10000   # KMC histogram x-axis cap; high-copy repeats pile up beyond this

mkdir -p "$OUTDIR"

echo '=== Step 1: k-mer profiling (the ceiling check, run BEFORE assembly) ==='
printf '%s\n%s\n' "$R1" "$R2" > "${OUTDIR}/reads.lst"
kmc -k${KMER} -t${THREADS} -m${MEM_GB} -ci1 -cs${MAX_KMER_COUNT} @"${OUTDIR}/reads.lst" "${OUTDIR}/kmcdb" "${OUTDIR}/"
kmc_tools transform "${OUTDIR}/kmcdb" histogram "${OUTDIR}/kmer.hist" -cx${MAX_KMER_COUNT}
genomescope2 -i "${OUTDIR}/kmer.hist" -o "${OUTDIR}/genomescope" -k ${KMER} -p ${PLOIDY}
echo "GenomeScope2 estimate written to ${OUTDIR}/genomescope/summary.txt"
echo "Inspect it: assembly total should land near this size; a het peak predicts fragmentation."

echo '=== Step 2: SPAdes isolate assembly (multi-k auto-selected; do NOT hand-pick k) ==='
# --isolate is tuned for high, even isolate coverage and does NOT run MismatchCorrector.
# Do not add --careful here: it is mutually exclusive with --isolate (SPAdes aborts if both are given).
spades.py --isolate -t ${THREADS} -m ${MEM_GB} -1 "$R1" -2 "$R2" -o "${OUTDIR}/spades"

ASM="${OUTDIR}/spades/contigs.fasta"   # report CONTIG metrics, not just scaffolds (scaffold gaps are estimated N-runs)
[ -f "$ASM" ] || { echo 'ERROR: assembly failed'; exit 1; }

echo '=== Step 3: contig-count sanity (biology, not failure) ==='
N_CONTIGS=$(grep -c '^>' "$ASM")
echo "Contigs: ${N_CONTIGS}"
echo "A clean isolate fragmenting into ~30-100 contigs is rRNA-operon / IS-element repeat structure,"
echo "not a tuning problem. For a single finished contig, use Unicycler hybrid or long reads."
echo "Next: NG50 + auN + BUSCO + Merqury QV (assembly-qc); compare total size to the GenomeScope estimate."
