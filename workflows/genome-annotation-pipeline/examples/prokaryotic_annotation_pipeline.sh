#!/bin/bash
# Reference: BRAKER3 3.0+, BUSCO 5.5+, Bakta 1.9+, Infernal 1.1+, InterProScan 5.66+, Prokka 1.14+, RepeatMasker 4.1+, RepeatModeler 2.0.4+, eggNOG-mapper 2.1+, pandas 2.2+, tRNAscan-SE 2.0+ | Verify API if version differs
# Complete prokaryotic genome annotation pipeline with Bakta
set -e

GENOME="assembly.fasta"
BAKTA_DB="/path/to/bakta_db"
PREFIX="my_genome"
LOCUS_TAG="MYORG"
THREADS=8
OUTDIR="bakta_annotation"

mkdir -p qc_reports    # NOT $OUTDIR: bakta exits if its --output dir already exists (unless --force)

echo "Step 1: Assembly QC"
# QUAST: basic assembly statistics
quast $GENOME -o qc_reports/quast --threads $THREADS

# BUSCO: genome completeness (bacteria_odb10 for bacteria)
busco -i $GENOME -l bacteria_odb10 -o busco_assembly --out_path qc_reports -m genome --cpu $THREADS

# CheckM2: the NON-NEGOTIABLE contamination gate before prokaryotic annotation. BUSCO does not
# measure contamination; >5% contam mixes two organisms' genes into one chimeric annotation.
checkm2 predict --input $GENOME --output-directory qc_reports/checkm2 --threads $THREADS

echo "CHECK: CheckM2 contamination < 5% AND completeness > 90% (and BUSCO/N50) before proceeding"

echo "Step 2: Bakta annotation"
bakta \
    --db $BAKTA_DB \
    --output $OUTDIR \
    --prefix $PREFIX \
    --locus-tag $LOCUS_TAG \
    --gram - \
    --translation-table 11 \
    --threads $THREADS \
    $GENOME
# --translation-table 11 for most bacteria; 4 for Mycoplasma/Spiroplasma (UGA = Trp) from
# the GTDB-Tk classification. Add --complete ONLY for finished replicons, not draft contigs.

echo "Step 3: Annotation QC"
# Count features from GFF3
echo "Feature counts:"
awk -F'\t' '$3 != "" && $1 !~ /^#/ {print $3}' $OUTDIR/${PREFIX}.gff3 | sort | uniq -c | sort -rn

# Check CDS count is reasonable for genome size
GENOME_SIZE=$(awk '/^[^>]/{total += length($0)} END{print total}' $GENOME)
CDS_COUNT=$(grep -c 'CDS' $OUTDIR/${PREFIX}.gff3 || echo 0)
TRNA_COUNT=$(grep -c 'tRNA' $OUTDIR/${PREFIX}.gff3 || echo 0)

echo ""
echo "Genome size: $GENOME_SIZE bp"
echo "CDS count: $CDS_COUNT"
echo "tRNA count: $TRNA_COUNT"

# Rule of thumb: ~1 gene per kbp for prokaryotes
EXPECTED_LOW=$((GENOME_SIZE / 1500))
EXPECTED_HIGH=$((GENOME_SIZE / 800))
echo "Expected CDS range: $EXPECTED_LOW - $EXPECTED_HIGH"

if [ "$CDS_COUNT" -lt "$EXPECTED_LOW" ] || [ "$CDS_COUNT" -gt "$EXPECTED_HIGH" ]; then
    echo "WARNING: CDS count outside expected range"
fi

if [ "$TRNA_COUNT" -lt 20 ]; then
    echo "WARNING: Fewer than 20 tRNAs detected"
fi

echo "Step 4: BUSCO on predicted proteins"
busco -i $OUTDIR/${PREFIX}.faa -l bacteria_odb10 -o busco_proteins --out_path qc_reports -m proteins --cpu $THREADS

echo ""
echo "Pipeline complete!"
echo "Outputs:"
echo "  GFF3: $OUTDIR/${PREFIX}.gff3"
echo "  GenBank: $OUTDIR/${PREFIX}.gbff"
echo "  Proteins: $OUTDIR/${PREFIX}.faa"
echo "  Nucleotides: $OUTDIR/${PREFIX}.ffn"
echo "  TSV summary: $OUTDIR/${PREFIX}.tsv"
