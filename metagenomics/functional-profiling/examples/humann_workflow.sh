#!/bin/bash
# Reference: HUMAnN 3.6+, MetaPhlAn 4.1+ | Verify API if version differs
# HUMAnN3 functional profiling: host-deplete and trim first (not shown), reuse the MetaPhlAn
# profile, normalize per sample, and KEEP UNMAPPED/UNINTEGRATED - they are the denominator.
set -euo pipefail

THREADS=8
OUTDIR="humann_results"
mkdir -p "$OUTDIR"

for fq in *.fastq.gz; do
    sample=$(basename "$fq" .fastq.gz)
    echo "Processing $sample..."
    # Do NOT pass --remove-temp-output: it deletes the MetaPhlAn profile you may want to reuse.
    humann --input "$fq" \
           --output "${OUTDIR}/${sample}" \
           --threads "$THREADS"
done

echo "Joining tables..."
humann_join_tables -i "$OUTDIR" -o "${OUTDIR}/merged_genefamilies.tsv" --file_name genefamilies
humann_join_tables -i "$OUTDIR" -o "${OUTDIR}/merged_pathabundance.tsv" --file_name pathabundance

# Normalize PER SAMPLE before cross-sample stats (RPK is depth-dependent); cpm preferred for models.
echo "Normalizing to CPM..."
humann_renorm_table -i "${OUTDIR}/merged_pathabundance.tsv" \
                    -o "${OUTDIR}/pathabundance_cpm.tsv" -u cpm

# Regroup to KEGG Orthologs; this adds an UNGROUPED row (analogue of UNINTEGRATED) - keep it.
echo "Regrouping to KEGG..."
humann_regroup_table -i "${OUTDIR}/merged_genefamilies.tsv" \
                     -g uniref90_ko \
                     -o "${OUTDIR}/merged_ko.tsv"

# Split into stratified (per-species) and unstratified (community totals incl. UNMAPPED/UNINTEGRATED).
# Run differential abundance on the UNSTRATIFIED file with a compositional method (MaAsLin2/ANCOM-BC).
echo "Splitting stratified tables..."
humann_split_stratified_table -i "${OUTDIR}/pathabundance_cpm.tsv" -o "$OUTDIR"
echo "Done. Results in ${OUTDIR}/ - remember these are functional POTENTIAL, not activity."
