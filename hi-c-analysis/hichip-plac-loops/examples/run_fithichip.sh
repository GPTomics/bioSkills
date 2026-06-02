#!/bin/bash
# Reference: FitHiChIP 11.0+, HiC-Pro 3.1+ | Verify API if version differs
# Calls FDR-controlled loops from a HiChIP/PLAC-seq library by writing a FitHiChIP
# config and running it. The caller, background, bias model, and FDR are set in the
# config (key=value), NOT on the command line - FitHiChIP_HiCPro.sh -C <config> is all.
set -euo pipefail

VALID_PAIRS=${1:-sample.allValidPairs.gz}   # HiC-Pro valid pairs for the HiChIP library
PEAK_FILE=${2:-chipseq_peaks.bed}           # anchors; prefer an INDEPENDENT ChIP-seq peak set
CHROM_SIZES=${3:-hg38.chrom.sizes}
OUT_DIR=${4:-fithichip_out}
PREFIX=${5:-sample}

# Mark sharpness drives the foreground/background pair:
#   broad mark (H3K27ac): IntType=3 (peak-to-all) + UseP2PBackgrnd=0 (loose)  -> sensitivity
#   sharp factor (CTCF):  IntType=1 (peak-to-peak) + UseP2PBackgrnd=1 (stringent) -> specificity
INT_TYPE=${INT_TYPE:-3}
USE_P2P_BG=${USE_P2P_BG:-0}

BIN_SIZE=5000        # 5kb: standard HiChIP anchor resolution (hichipper effective ~2.5kb)
LOW_DIST=20000       # 20kb floor: below this, self-ligation/diagonal dominate, not loops
UPP_DIST=2000000     # 2Mb ceiling: loops beyond this are rare and noise-dominated at HiChIP depth
BIAS_TYPE=1          # 1=coverage-bias regression (default); 2=ICE-bias regression
QVALUE=0.01          # FDR cutoff for significant loops

CONFIG=$(mktemp config_fithichip.XXXX)
cat > "$CONFIG" <<EOF
ValidPairs=${VALID_PAIRS}
PeakFile=${PEAK_FILE}
ChrSizeFile=${CHROM_SIZES}
OutDir=${OUT_DIR}/
PREFIX=${PREFIX}
BINSIZE=${BIN_SIZE}
LowDistThr=${LOW_DIST}
UppDistThr=${UPP_DIST}
IntType=${INT_TYPE}
UseP2PBackgrnd=${USE_P2P_BG}
BiasType=${BIAS_TYPE}
MergeInt=1
QVALUE=${QVALUE}
EOF

bash FitHiChIP_HiCPro.sh -C "$CONFIG"

# Significant loops: ${OUT_DIR}/${PREFIX}.interactions_FitHiC_Q${QVALUE}.bed
# (and ..._MergeNearContacts.bed because MergeInt=1)
echo "Loops in ${OUT_DIR}/${PREFIX}.interactions_FitHiC_Q${QVALUE}.bed"
rm -f "$CONFIG"
