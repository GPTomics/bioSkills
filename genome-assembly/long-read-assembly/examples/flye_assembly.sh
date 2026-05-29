#!/bin/bash
# Reference: Flye 2.9+, Porechop_ABI 0.5+, minimap2 2.26+ | Verify API if version differs
# Era-matched noisy-long-read assembly with Flye: de-chimerize, pick the flag from the
# basecaller era, assemble, then flag the haplotig-duplication red flag and the polishing handoff.
# Usage: flye_assembly.sh <reads.fq.gz> <outdir> <era> [expected_genome_size] [threads]
#   era = r10 | r9sup | r9raw | clr   (the basecaller era IS an assembly parameter, not metadata)
set -euo pipefail

READS=$1
OUTDIR=$2
ERA=$3
SIZE=${4:-}            # expected genome size (e.g. 4.6m, 3g); optional in recent Flye, used for the dup check
THREADS=${5:-16}

# The most dangerous mistake is telling the assembler the reads are noisier than they are:
# --nano-raw on R10/Dorado data over-collapses repeats and RAISES N50 as the assembly gets worse.
case "$ERA" in
  r10)   MODE='--nano-hq' ;;                         # ONT R10.4.1 / Dorado HAC-SUP (Q20+)
  r9sup) MODE='--nano-hq --read-error 0.05' ;;       # ONT R9.4.1 SUP: hq mode with the error floor raised
  r9raw) MODE='--nano-raw' ;;                         # legacy ONT R9 fast/HAC only
  clr)   MODE='--pacbio-raw' ;;                       # PacBio CLR (legacy; polish hard afterward)
  *)     echo "Unknown era '$ERA' (use r10|r9sup|r9raw|clr)"; exit 1 ;;
esac

echo "=== Noisy long-read assembly (era=$ERA, mode=$MODE) ==="

# 1. De-chimerize internal adapters (ONT-specific): an untrimmed internal adapter becomes a
#    structural mis-join that polishing cannot fix. ab initio detection; splits chimeric reads.
TRIMMED=${OUTDIR}/trimmed.fq.gz
mkdir -p "$OUTDIR"
if [ "$ERA" != "clr" ]; then
    porechop_abi -abi -i "$READS" -o "$TRIMMED" -t "$THREADS"
else
    TRIMMED=$READS    # PacBio removes adapters/scraps upstream on the instrument
fi

# 2. Assemble. --genome-size is optional in recent Flye but required with --asm-coverage.
SIZE_ARG=""
[ -n "$SIZE" ] && SIZE_ARG="--genome-size $SIZE"
flye $MODE "$TRIMMED" $SIZE_ARG --out-dir "${OUTDIR}/flye" --threads "$THREADS"

ASM=${OUTDIR}/flye/assembly.fasta
echo ""
echo "=== Per-contig stats (length / coverage / circular) ==="
cat "${OUTDIR}/flye/assembly_info.txt"

# 3. Haplotig false-duplication red flag: assembled size 1.5-2x the expected genome size means the
#    assembler kept both haplotypes as primary contigs (confirm with a half-coverage depth peak and
#    inflated BUSCO-Duplicated, then run purge_dups). 1.5 is the conventional red-flag floor.
DUP_FLAG=1.5
if [ -n "$SIZE" ]; then
    EXP_BP=$(python3 -c "s='$SIZE'.lower(); u={'k':1e3,'m':1e6,'g':1e9}; print(int(float(s[:-1])*u[s[-1]]) if s[-1] in u else int(s))")
    ASM_BP=$(grep -v '^>' "$ASM" | tr -d '\n' | wc -c | tr -d ' ')
    RATIO=$(python3 -c "print(f'{$ASM_BP/$EXP_BP:.2f}')")
    echo ""
    echo "Assembled ${ASM_BP} bp vs expected ${EXP_BP} bp (ratio ${RATIO})"
    python3 -c "import sys; sys.exit(0 if $RATIO < $DUP_FLAG else 1)" \
        || echo "RED FLAG: ratio >= ${DUP_FLAG} suggests haplotig duplication -> run purge_dups (inspect coverage cutoffs first)"
fi

echo ""
echo "Assembly: $ASM"
echo "NOT FINISHED: raw long-read consensus is low-QV (indels frameshift genes)."
echo "Next: polish (assembly-polishing) then measure QV with Merqury (assembly-qc)."
