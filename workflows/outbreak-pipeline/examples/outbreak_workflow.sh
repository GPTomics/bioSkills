#!/bin/bash
# Outbreak investigation pipeline: MLST typing, AMR, phylodynamics, transmission
# Requires: mlst, abricate, snippy, iqtree2, treetime
set -e

ISOLATES_DIR="isolates"
REFERENCE="reference.gbk"
METADATA="metadata.tsv"  # name\tdate format
OUTDIR="outbreak_results"
THREADS=8

# Generation time parameters (adjust for pathogen)
# Bacteria: ~7-14 days; Viruses: 4-7 days; TB: ~1 year
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

# Count SNPs
n_snps=$(grep -v "^>" ${OUTDIR}/alignment/core.aln | tr -d '\n-' | wc -c)
n_sites=$(grep -v "^>" ${OUTDIR}/alignment/core.aln | head -1 | tr -d '\n' | wc -c)
echo "Core alignment: ${n_snps} variable sites in ${n_sites} bp"

# === Step 3: Phylogenetic Tree ===
echo ""
echo "=== Step 3: Phylogenetic Tree ==="
iqtree2 -s ${OUTDIR}/alignment/core.aln \
        -m GTR+G \
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
    --aln ${OUTDIR}/alignment/core.aln \
    --dates ${METADATA} \
    --outdir ${OUTDIR}/phylo/treetime_output \
    --coalescent skyline \
    --clock-filter 3 \
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

# Run inference
res <- inferTTree(tree, dateT = dateT,
                   w.shape = w_shape, w.scale = w_scale,
                   mcmcIterations = 10000,
                   startNeg = 1, startPi = 0.5)

# Save results
saveRDS(res, "${OUTDIR}/transmission/transphylo_result.rds")

# Extract median tree
medTree <- medTTree(res)

# Plot transmission tree
pdf("${OUTDIR}/transmission/transmission_tree.pdf", width=12, height=10)
plotTTree(medTree)
dev.off()

# Who infected whom matrix
wiw <- computeMatWIW(res)
write.csv(wiw, "${OUTDIR}/transmission/who_infected_whom.csv")

# R0 estimate
offspring <- getOffspringMulti(res)
cat("\n=== R0 Estimate ===\n")
cat("Mean R0:", mean(offspring), "\n")
cat("95% CI:", quantile(offspring, 0.025), "-", quantile(offspring, 0.975), "\n")

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
