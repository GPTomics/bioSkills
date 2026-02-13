#!/bin/bash
# Reference: BRAKER3 3.0+, BUSCO 5.5+, Bakta 1.9+, Infernal 1.1+, InterProScan 5.66+, Prokka 1.14+, RepeatMasker 4.1+, RepeatModeler 2.0+, eggNOG-mapper 2.1+, pandas 2.2+, tRNAscan-SE 2.0+ | Verify API if version differs
# Complete eukaryotic genome annotation pipeline
# RepeatModeler/RepeatMasker -> BRAKER3 -> eggNOG-mapper -> Infernal/tRNAscan-SE
set -e

GENOME="assembly.fasta"
RNASEQ_BAM="rnaseq_sorted.bam"
PROTEINS="orthodb_proteins.fa"
BUSCO_LINEAGE="eukaryota_odb10"
THREADS=8

mkdir -p qc_reports repeat_out braker_out functional_out ncrna_out

echo "Step 0: Assembly QC"
quast $GENOME -o qc_reports/quast --threads $THREADS --eukaryote
busco -i $GENOME -l $BUSCO_LINEAGE -o qc_reports/busco_assembly -m genome --cpu $THREADS
echo "CHECK: Verify assembly quality before proceeding"

echo "Step 1: Repeat modeling and masking"
# Build species-specific repeat library (~hours to days depending on genome size)
RepeatModeler -database mygenome -pa $THREADS -LTRStruct

# Mask the genome with soft-masking (-xsmall) for downstream gene prediction
RepeatMasker \
    -lib mygenome-families.fa \
    -pa $THREADS \
    -xsmall \
    -gff \
    -dir repeat_out \
    $GENOME

echo "Repeat masking summary:"
tail -20 repeat_out/$(basename $GENOME).tbl

# QC: Check repeat content
TOTAL_MASKED=$(grep 'total interspersed' repeat_out/$(basename $GENOME).tbl | awk '{print $NF}')
echo "Total interspersed repeat content: $TOTAL_MASKED"
echo "CHECK: Verify repeat content is within expected range for taxon"
echo "  Vertebrate: 30-60%, Insect: 15-45%, Plant: 20-85%, Fungus: 3-20%"

MASKED_GENOME="repeat_out/$(basename $GENOME).masked"

echo "Step 2: Gene prediction with BRAKER3"
# BRAKER3 combines GeneMark-ETP, AUGUSTUS, and TSEBRA
# Uses both RNA-seq and protein evidence for best results
braker.pl \
    --genome=$MASKED_GENOME \
    --bam=$RNASEQ_BAM \
    --prot_seq=$PROTEINS \
    --softmasking \
    --threads $THREADS \
    --species=my_species \
    --gff3 \
    --workingdir=braker_out

# Count predicted genes
GENE_COUNT=$(grep -c 'gene' braker_out/braker.gff3 | head -1)
EXON_COUNT=$(grep -c 'exon' braker_out/braker.gff3 | head -1)
echo "Predicted genes: $GENE_COUNT"
echo "Predicted exons: $EXON_COUNT"

echo "Step 2b: BUSCO QC on predicted proteins"
busco -i braker_out/braker.aa -l $BUSCO_LINEAGE -o qc_reports/busco_proteins -m protein --cpu $THREADS

echo "CHECK: BUSCO completeness should be > 90%"

echo "Step 3: Functional annotation with eggNOG-mapper"
emapper.py \
    -i braker_out/braker.aa \
    --output functional_out/eggnog \
    --cpu $THREADS \
    -m diamond \
    --tax_scope auto \
    --go_evidence non-electronic \
    --target_orthologs all \
    --seed_ortholog_evalue 1e-5 \
    --override

# Count annotated genes
TOTAL_PROTEINS=$(grep -c '>' braker_out/braker.aa)
ANNOTATED=$(grep -v '^#' functional_out/eggnog.emapper.annotations | grep -v '^$' | wc -l)
PCT_ANNOTATED=$((ANNOTATED * 100 / TOTAL_PROTEINS))
echo "Functional annotation: $ANNOTATED / $TOTAL_PROTEINS ($PCT_ANNOTATED%)"

if [ "$PCT_ANNOTATED" -lt 60 ]; then
    echo "WARNING: < 60% annotated. Try broader --tax_scope or add InterProScan"
fi

echo "Step 3b: InterProScan for domain annotation"
interproscan.sh \
    -i braker_out/braker.aa \
    -o functional_out/interpro_results.tsv \
    -f tsv,gff3 \
    -goterms \
    -pa \
    -cpu $THREADS

echo "Step 4: ncRNA annotation"
# tRNAscan-SE for tRNA genes (-E for eukaryotic)
tRNAscan-SE \
    -E \
    --thread $THREADS \
    -o ncrna_out/trna_results.txt \
    --gff ncrna_out/trna.gff \
    $GENOME

# Infernal for Rfam-based ncRNA annotation
cmscan \
    --cpu $THREADS \
    --tblout ncrna_out/rfam_results.tbl \
    --fmt 2 \
    --clanin Rfam.clanin \
    Rfam.cm \
    $GENOME

TRNA_COUNT=$(grep -c 'tRNA' ncrna_out/trna.gff || echo 0)
echo "tRNAs detected: $TRNA_COUNT"

echo "Step 5: Merge and validate annotations"
# AGAT for GFF merging (if available)
if command -v agat_sp_merge_annotations.pl &> /dev/null; then
    agat_sp_merge_annotations.pl \
        --gff braker_out/braker.gff3 \
        --gff ncrna_out/trna.gff \
        -o final_annotation.gff3

    agat_sp_statistics.pl \
        --gff final_annotation.gff3 \
        -o annotation_stats.txt
else
    echo "AGAT not available. Merge annotations manually or install AGAT."
    cp braker_out/braker.gff3 final_annotation.gff3
fi

echo ""
echo "Pipeline complete!"
echo "Outputs:"
echo "  Gene models: braker_out/braker.gff3"
echo "  Proteins: braker_out/braker.aa"
echo "  eggNOG: functional_out/eggnog.emapper.annotations"
echo "  InterPro: functional_out/interpro_results.tsv"
echo "  tRNAs: ncrna_out/trna.gff"
echo "  Rfam ncRNAs: ncrna_out/rfam_results.tbl"
echo "  Repeat masking: repeat_out/"
echo "  QC: qc_reports/"
