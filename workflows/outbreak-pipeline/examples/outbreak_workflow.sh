#!/bin/bash
# Reference: AMRFinderPlus 3.12+, BioPython 1.83+, IQ-TREE 2.2+, Nextclade 3.3+, TreeTime 0.11+, matplotlib 3.8+, mlst 2.23+, pandas 2.2+ | Verify API if version differs
# Outbreak investigation pipeline: MLST typing, AMR, phylodynamics, transmission
# Requires: mlst, abricate, snippy, iqtree2, treetime
set -e

ISOLATES_DIR="isolates"
REFERENCE="reference.gbk"
METADATA="metadata.tsv"  # name\tdate format
OUTDIR="outbreak_results"
THREADS=8

# Generation-time Gamma SCALE parameter in days (mean = shape*scale; w_shape=2 below, so
# GEN_TIME_DAYS=7 gives a ~14-day mean generation time). Set the SCALE for your pathogen and
# cite the literature. Bacteria: ~7-14 day mean; Viruses: 4-7 day; TB: ~1 year.
GEN_TIME_DAYS=7

mkdir -p ${OUTDIR}/{mlst,amr,alignment,phylo,transmission,qc}

echo "=== Outbreak Investigation Pipeline ==="
echo "Isolates: ${ISOLATES_DIR}"
echo "Reference: ${REFERENCE}"
echo "Metadata: ${METADATA}"

# === Step 1a: MLST Typing ===
echo ""
echo "=== Step 1a: MLST Typing ==="
for fasta in ${ISOLATES_DIR}/*.fasta; do
    sample=$(basename $fasta .fasta)
    echo "MLST: ${sample}"
    mlst $fasta >> ${OUTDIR}/mlst/all_mlst.tsv 2>/dev/null
done

# Check clonality
echo ""
echo "MLST Results:"
cut -f2,3 ${OUTDIR}/mlst/all_mlst.tsv | sort | uniq -c | sort -rn
n_sts=$(cut -f3 ${OUTDIR}/mlst/all_mlst.tsv | sort -u | wc -l)
if [ $n_sts -gt 1 ]; then
    echo "WARNING: Multiple sequence types detected - outbreak may involve multiple clones"
fi

# === Step 1b: AMR Detection (parallel with MLST) ===
# NOTE: this demo uses abricate for brevity. The skill's preferred path is AMRFinderPlus --organism
# (species-aware, reports allele identity not just gene family) unified across tools with hAMRonization.
echo ""
echo "=== Step 1b: AMR Detection ==="
for fasta in ${ISOLATES_DIR}/*.fasta; do
    sample=$(basename $fasta .fasta)
    echo "AMR: ${sample}"
    abricate --db ncbi $fasta > ${OUTDIR}/amr/${sample}.ncbi.tsv 2>/dev/null
    abricate --db card $fasta > ${OUTDIR}/amr/${sample}.card.tsv 2>/dev/null
done

# Summary matrix
abricate --summary ${OUTDIR}/amr/*.ncbi.tsv > ${OUTDIR}/amr/amr_summary_ncbi.tsv
abricate --summary ${OUTDIR}/amr/*.card.tsv > ${OUTDIR}/amr/amr_summary_card.tsv
echo "AMR summaries: ${OUTDIR}/amr/"

# === Step 2: Core Genome Alignment ===
echo ""
echo "=== Step 2: Core Genome Alignment ==="
for fasta in ${ISOLATES_DIR}/*.fasta; do
    sample=$(basename $fasta .fasta)
    echo "Snippy: ${sample}"
    snippy --outdir ${OUTDIR}/alignment/snippy_${sample} \
           --ref ${REFERENCE} \
           --ctgs $fasta \
           --cpus ${THREADS} \
           --quiet
done

# Core SNP alignment
echo "Building core alignment..."
cd ${OUTDIR}/alignment
snippy-core --ref ../../${REFERENCE} snippy_*
cd ../..

# === Recombination masking (MANDATORY for recombining bacteria; skip ONLY for clonal Mtb) ===
# Gubbins needs core.full.aln (FULL positions incl. invariant), NOT core.aln (variable-only), to
# estimate background SNP density. Skipping this inflates the clock rate 2-5x and fabricates links.
echo "Masking recombination (Gubbins)..."
run_gubbins.py --prefix ${OUTDIR}/alignment/gubbins ${OUTDIR}/alignment/core.full.aln
# Gubbins emits the recombination-filtered polymorphic-sites alignment used for the tree:
MASKED_ALN=${OUTDIR}/alignment/gubbins.filtered_polymorphic_sites.fasta

# === Step 3: Phylogenetic Tree (on the masked alignment; +ASC because input is variable-sites-only) ===
echo ""
echo "=== Step 3: Phylogenetic Tree ==="
iqtree2 -s ${MASKED_ALN} \
        -m GTR+G+ASC \
        -bb 1000 \
        -nt AUTO \
        --prefix ${OUTDIR}/phylo/outbreak \
        -quiet

echo "Tree: ${OUTDIR}/phylo/outbreak.treefile"

# === Step 4: Phylodynamics with TreeTime ===
echo ""
echo "=== Step 4: Phylodynamics ==="

# Check metadata exists
if [ ! -f "${METADATA}" ]; then
    echo "ERROR: Metadata file ${METADATA} not found"
    echo "Create a TSV with columns: name, date (YYYY-MM-DD)"
    exit 1
fi

treetime \
    --tree ${OUTDIR}/phylo/outbreak.treefile \
    --aln ${MASKED_ALN} \
    --dates ${METADATA} \
    --outdir ${OUTDIR}/phylo/treetime_output \
    --coalescent skyline \
    --clock-filter 4 \
    --confidence

# Check temporal signal
if [ -f "${OUTDIR}/phylo/treetime_output/root_to_tip_regression.pdf" ]; then
    echo "Temporal signal plot: ${OUTDIR}/phylo/treetime_output/root_to_tip_regression.pdf"
fi

# Extract clock rate
if [ -f "${OUTDIR}/phylo/treetime_output/dates.tsv" ]; then
    echo ""
    echo "Dated tree: ${OUTDIR}/phylo/treetime_output/timetree.nexus"
fi

# === Step 5: Transmission Inference ===
echo ""
echo "=== Step 5: Transmission Inference ==="

# Get latest date from metadata for dateT
latest_date=$(cut -f2 ${METADATA} | tail -n +2 | sort | tail -1)
# Convert to decimal year (approximate)
year=$(echo $latest_date | cut -d'-' -f1)
month=$(echo $latest_date | cut -d'-' -f2)
day=$(echo $latest_date | cut -d'-' -f3)
decimal_date=$(echo "scale=4; $year + ($month - 1) / 12 + $day / 365" | bc)

# Generation time in years
gen_time_years=$(echo "scale=6; ${GEN_TIME_DAYS} / 365" | bc)

cat > ${OUTDIR}/transmission/run_transphylo.R << EOF
library(TransPhylo)
library(ape)

# Load dated tree
tree <- read.nexus("${OUTDIR}/phylo/treetime_output/timetree.nexus")

# Parameters
dateT <- ${decimal_date}
w_shape <- 2
w_scale <- ${gen_time_years}

cat("Running TransPhylo...\n")
cat("Date T:", dateT, "\n")
cat("Generation time: shape=", w_shape, ", scale=", w_scale, " years\n")

# TransPhylo needs a ptree (dated phylo + last-sample date), not a raw ape phylo.
ptree <- ptreeFromPhylo(tree, dateLastSample = dateT)

# Run inference
res <- inferTTree(ptree, dateT = dateT,
                   w.shape = w_shape, w.scale = w_scale,
                   mcmcIterations = 10000,
                   startNeg = 1, startPi = 0.5)

# Save results
saveRDS(res, "${OUTDIR}/transmission/transphylo_result.rds")

# Extract median coloured transmission tree (ctree)
med_ctree <- medTTree(res)

# Plot transmission tree
pdf("${OUTDIR}/transmission/transmission_tree.pdf", width=12, height=10)
plotCTree(med_ctree)
dev.off()

# Who infected whom matrix (same 0.5 burn-in as the R_e estimate below)
wiw <- computeMatWIW(res, burnin = 0.5)
write.csv(wiw, "${OUTDIR}/transmission/who_infected_whom.csv")

# Effective reproduction number (R_e): getOffspringDist is per-case, average across sampled hosts
# (names come from the ptree; res has no \$ttree\$nam field)
offspring <- sapply(ptree\$nam, function(k) mean(getOffspringDist(res, k = k, burnin = 0.5)))
cat("\n=== R_e Estimate ===\n")
cat("Mean R_e:", mean(offspring), "\n")
# Across-host 2.5-97.5% spread of per-host mean offspring, NOT a posterior credible interval.
cat("Across-host 2.5-97.5%:", quantile(offspring, 0.025), "-", quantile(offspring, 0.975), "\n")

# Summary
cat("\n=== Transmission Summary ===\n")
cat("Number of sampled cases:", length(tree\$tip.label), "\n")
EOF

Rscript ${OUTDIR}/transmission/run_transphylo.R

# === Summary ===
echo ""
echo "=== Pipeline Complete ==="
echo ""
echo "Results:"
echo "  MLST types:          ${OUTDIR}/mlst/all_mlst.tsv"
echo "  AMR summary:         ${OUTDIR}/amr/amr_summary_ncbi.tsv"
echo "  Core alignment:      ${OUTDIR}/alignment/core.aln"
echo "  ML tree:             ${OUTDIR}/phylo/outbreak.treefile"
echo "  Dated tree:          ${OUTDIR}/phylo/treetime_output/timetree.nexus"
echo "  Transmission tree:   ${OUTDIR}/transmission/transmission_tree.pdf"
echo "  Transmission matrix: ${OUTDIR}/transmission/who_infected_whom.csv"
