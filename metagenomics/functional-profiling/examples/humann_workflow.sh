#!/bin/bash
# HUMAnN3 functional profiling workflow

THREADS=8
OUTDIR="humann_results"
mkdir -p $OUTDIR

for fq in *.fastq.gz; do
    sample=$(basename $fq .fastq.gz)
    echo "Processing $sample..."

    humann --input $fq \
           --output ${OUTDIR}/${sample} \
           --threads $THREADS \
           --remove-temp-output
done

echo "Joining tables..."
humann_join_tables -i $OUTDIR -o ${OUTDIR}/merged_genefamilies.tsv --file_name genefamilies
humann_join_tables -i $OUTDIR -o ${OUTDIR}/merged_pathabundance.tsv --file_name pathabundance

echo "Normalizing..."
humann_renorm_table -i ${OUTDIR}/merged_pathabundance.tsv \
                    -o ${OUTDIR}/pathabundance_relab.tsv -u relab

echo "Regrouping to KEGG..."
humann_regroup_table -i ${OUTDIR}/merged_genefamilies.tsv \
                     -g uniref90_ko \
                     -o ${OUTDIR}/merged_ko.tsv

echo "Splitting stratified tables..."
humann_split_stratified_table -i ${OUTDIR}/pathabundance_relab.tsv -o ${OUTDIR}

echo "Done! Results in ${OUTDIR}/"
