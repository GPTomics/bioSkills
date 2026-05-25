#!/bin/bash
# Reference: NCBI BLAST+ 2.15+, HMMER 3.4+, MMseqs2 15+ | Verify API if version differs
# Compare three iterative profile search strategies: psiblast (3 iter), jackhmmer (3 iter), MMseqs2 (3 iter at -s 7.5).

set -euo pipefail

QUERY="${1:-distant_protein.fa}"
DB="${2:-uniref90}"   # blast/HMMER FASTA path; MMseqs2 will build its own DB

ITERATIONS=3
INCL_E=0.002          # Tighter than default 0.005 -- reduces paralog drift
EVALUE=1e-5
THREADS=8

echo "=== PSI-BLAST ${ITERATIONS} iterations, -inclusion_ethresh ${INCL_E} ==="
psiblast -query "${QUERY}" -db "${DB}" \
         -num_iterations "${ITERATIONS}" \
         -inclusion_ethresh "${INCL_E}" \
         -evalue "${EVALUE}" \
         -num_threads "${THREADS}" \
         -out_pssm psiblast.pssm.asn \
         -out_ascii_pssm psiblast.pssm.txt \
         -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore" \
         -out psiblast.tsv
psiblast_hits=$(wc -l < psiblast.tsv)
echo "  psiblast: ${psiblast_hits} HSPs"

echo
echo "=== jackhmmer ${ITERATIONS} iterations ==="
jackhmmer -N "${ITERATIONS}" \
          --tblout jackhmmer.tbl \
          --chkhmm jackhmmer.hmm \
          --cpu "${THREADS}" \
          "${QUERY}" "${DB}" > jackhmmer.txt
jackhmmer_hits=$(grep -cv '^#' jackhmmer.tbl)
echo "  jackhmmer: ${jackhmmer_hits} hits"

echo
echo "=== MMseqs2 ${ITERATIONS} iterations at -s 7.5 ==="
mkdir -p mmseqs_tmp
mmseqs createdb "${QUERY}" mmseqs_queryDB
mmseqs createdb "${DB}" mmseqs_targetDB
mmseqs search mmseqs_queryDB mmseqs_targetDB mmseqs_resultDB mmseqs_tmp \
       --num-iterations "${ITERATIONS}" -s 7.5 -e "${EVALUE}" --threads "${THREADS}"
mmseqs convertalis mmseqs_queryDB mmseqs_targetDB mmseqs_resultDB mmseqs.m8 \
       --format-output query,target,fident,alnlen,evalue,bits
mmseqs_hits=$(wc -l < mmseqs.m8)
echo "  mmseqs2: ${mmseqs_hits} hits"

echo
echo "=== Summary ==="
printf 'psiblast :  %6d HSPs (slowest)\n' "${psiblast_hits}"
printf 'jackhmmer:  %6d hits  (sensitivity often best)\n' "${jackhmmer_hits}"
printf 'mmseqs2  :  %6d hits  (fastest at comparable sensitivity)\n' "${mmseqs_hits}"
echo
echo "Caveat: counts are not directly comparable -- different scoring + reporting."
echo "Cross-validate top hits across methods; trust hits that appear in 2+ methods."
