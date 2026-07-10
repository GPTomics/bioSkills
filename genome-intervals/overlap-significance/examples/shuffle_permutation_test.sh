#!/bin/bash
# Reference: bedtools 2.31+ | Verify API if version differs
# Permutation null for interval-overlap significance with bedtools shuffle + jaccard.
# Triages with fisher first, then builds an empirical null restricted to an accessible
# workspace and outside the ENCODE blacklist - the matched null, not a uniform-random one.
set -euo pipefail

QUERY=peaks.bed
ANNOTATION=enhancers.bed
GENOME=genome.txt          # chrom<TAB>size; defines chromosome order and ends
WORKSPACE=accessible.bed   # the universe: regions the query could have come from (NOT the whole genome)
BLACKLIST=blacklist.bed    # ENCODE blacklist + assembly gaps; excluded from placement
N_PERMUTATIONS=1000        # >=1000 gives a stable empirical p down to ~0.001

sort -k1,1 -k2,2n "$QUERY" > query.sorted.bed
sort -k1,1 -k2,2n "$ANNOTATION" > annotation.sorted.bed

echo "== fisher triage (weak analytic null - screen only) =="
bedtools fisher -a query.sorted.bed -b annotation.sorted.bed -g "$GENOME"

observed=$(bedtools jaccard -a query.sorted.bed -b annotation.sorted.bed | awk 'NR==2{print $3}')
echo "observed jaccard: $observed"

: > null_jaccards.txt
for _ in $(seq 1 "$N_PERMUTATIONS"); do
    bedtools shuffle -i query.sorted.bed -g "$GENOME" -incl "$WORKSPACE" -excl "$BLACKLIST" -chrom \
        | sort -k1,1 -k2,2n \
        | bedtools jaccard -a - -b annotation.sorted.bed \
        | awk 'NR==2{print $3}' >> null_jaccards.txt
done

# Empirical one-sided p with the +1 correction (Phipson & Smyth 2010): a permutation p is never 0.
awk -v obs="$observed" -v n="$N_PERMUTATIONS" '
    $1 >= obs {hits++}
    END {printf "empirical p = %.4g  (hits=%d / N=%d)\n", (hits + 1) / (n + 1), hits, n}' null_jaccards.txt
