#!/bin/bash
# Download reference panels for imputation

set -euo pipefail

PANEL=${1:-1000G}
BUILD=${2:-GRCh38}
OUTPUT_DIR=${3:-reference_panels}

mkdir -p $OUTPUT_DIR

echo "=== Downloading Reference Panel: $PANEL ==="

case $PANEL in
    "1000G")
        echo "Downloading 1000 Genomes Phase 3..."
        # Download from IGSR
        for chr in {1..22}; do
            wget -P $OUTPUT_DIR \
                "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz"
        done
        ;;

    "HRC")
        echo "HRC panel requires application - see http://www.haplotype-reference-consortium.org/"
        ;;

    "TOPMed")
        echo "TOPMed requires dbGaP access - see https://topmed.nhlbi.nih.gov/"
        ;;

    *)
        echo "Unknown panel: $PANEL"
        echo "Available: 1000G, HRC, TOPMed"
        exit 1
        ;;
esac

echo "Download complete. Files in: $OUTPUT_DIR"
