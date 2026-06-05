#!/usr/bin/env bash
# Reference: NetMHCIIpan 4.3+, MixMHC2pred 2.0+ | Verify API if version differs
# Predict HLA class II presentation with NetMHCIIpan-4.3 (default EL %Rank, -BA adds
# affinity) and MixMHC2pred-2.0. Allele nomenclature differs by tool AND isotype:
# DR is single-chain; DQ/DP are alpha/beta heterodimers needing both chains.
set -euo pipefail

PEPTIDES=${1:-peptides.txt}      # one peptide per line, 12-25mers
ANTIGEN=${2:-antigen.fasta}

# DR (single-chain): the beta allele names the molecule.
netMHCIIpan -f "$PEPTIDES" -inptype 1 -a DRB1_0101,DRB1_1501 -BA -xls -xlsfile dr_out.tsv

# DQ / DP heterodimers: BOTH chains, hyphen-joined. Use documented cis pairings only;
# do NOT expand to every DQA1 x DQB1 combination (creates non-existent molecules).
netMHCIIpan -f "$ANTIGEN" -inptype 0 -length 15 \
    -a HLA-DQA10501-DQB10201,HLA-DPA10103-DPB10401 -xls -xlsfile dqdp_out.tsv

# MixMHC2pred: chain-underscore names, DOUBLE underscore between chains, space-separated.
MixMHC2pred -i "$PEPTIDES" -o mixmhc2pred_out.txt \
    -a DRB1_15_01 DRB5_01_01 DPA1_02_01__DPB1_01_01

# Class II %Rank cutoffs are LOOSER than class I: strong <= 1%, weak <= 5%.
echo "Strong class II binders (Rank <= 1%) in dr_out.tsv:"
awk -F'\t' 'NR>1 && $NF<=1.0 {print}' dr_out.tsv | head
