# Genome Distance and Species Delineation - Usage Guide

## Overview

Prokaryote species delineation has shifted from 16S rRNA identity (insufficient at < 98.7%) to whole-genome **ANI at 95%** (Jain 2018 Nat Commun 9:5114). The modern operational standard for taxonomic classification is **GTDB-Tk** (Chaumeil 2020/2022 Bioinformatics 38:5315), using GTDB taxonomy + ANI radius. **skani** (Shaw & Yu 2023 Nat Methods 20:1661) has replaced FastANI as the default ANI tool in GTDB-Tk 2.4+, being 20-30x faster while maintaining accuracy.

The 95% ANI threshold is robust but genus-specific (Parks 2018 Nat Biotech 36:996 demonstrated radius varies). For novel-species naming, **dDDH** (digital DNA-DNA hybridization) via TYGS / GGDC remains the gold standard.

## Prerequisites

```bash
# Modern fast ANI
conda install -c bioconda skani

# Standard FastANI
conda install -c bioconda fastani

# Multi-method
pip install pyani-plus

# Mash for fast clustering
conda install -c bioconda mash dashing2

# GTDB-Tk for taxonomy
conda install -c bioconda gtdbtk
# Download GTDB release (must match GTDB-Tk version)
wget https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_r220_data.tar.gz
tar xf gtdbtk_r220_data.tar.gz
export GTDBTK_DATA_PATH=$PWD/release220

# MAG QC
conda install -c bioconda checkm2
```

## Quick Start

Tell the AI agent what to do:
- "Compute pairwise ANI between these 50 bacterial genomes using skani"
- "Classify these MAGs to GTDB taxonomy using GTDB-Tk"
- "Apply the 95% ANI + AF >= 0.5 species delineation threshold to identify same-species clusters"
- "Generate Mash distance matrix for fast clustering of 1000 genomes"

## Example Prompts

### Bacterial Strain Comparison

> "Compute ANI between these 30 E. coli strains using skani. Apply the 95% ANI species delineation threshold with alignment fraction >= 0.5. Classify each strain to GTDB-r220 taxonomy using GTDB-Tk. Cross-validate with FastANI for borderline pairs (94-96% ANI)."

### MAG Taxonomy

> "Classify these 200 metagenome-assembled genomes (MAGs) to GTDB taxonomy. First filter MAGs with CheckM2 (>= 70% completeness, < 5% contamination). Then run GTDB-Tk classify_wf. Report classification at each rank (kingdom, phylum, class, order, family, genus, species) with completeness metrics."

### Novel Species Delineation

> "I have a putative novel bacterial species. Compare it against the closest type strain using TYGS for dDDH calculation. Cross-validate with skani ANI (expecting < 95%) and verify all three signals (ANI, dDDH, AAI) indicate distinct species. Report the closest established type strain and dDDH percentage."

## What the Agent Will Do

1. **Quality-check inputs**: Compleasm completeness per genome; CheckM2 for MAGs
2. **Sketch genomes** for skani / FastANI (or compute on-the-fly)
3. **Compute ANI** with skani (fast) or FastANI (slower; legacy comparison)
4. **Apply species delineation**: ANI >= 95% + AF >= 0.5
5. **Cross-validate** with OrthoANI for borderline cases (94-96% ANI)
6. **Run GTDB-Tk** for taxonomic classification
7. **Report**: pairwise ANI matrix, GTDB taxonomy, alignment fractions
8. **For novel species**: TYGS dDDH via web; type strain comparison
9. **Caveats**: genus-specific ANI radius, contamination effects, ANI saturation < 75%

## Tips

- skani (Shaw 2023) is the modern default; 20-30x faster than FastANI; default in GTDB-Tk 2.4+
- For ANI matrix of 1000+ genomes, use `skani triangle` all-vs-all mode
- ANI saturates below 75%; switch to AAI for distant comparisons
- Alignment fraction (AF) >= 0.5 required for valid species call (Jain 2018 standard)
- GTDB-Tk uses skani; pin GTDB-Tk + GTDB release version compatibility
- For MAGs, CheckM2 pre-filter is essential (>= 70% completeness, < 5% contamination)
- dDDH via TYGS or GGDC web service is gold standard for novel species naming
- Mash distance is fast but approximate (0.05 ~ 95% ANI heuristic, not exact)
- For sub-species / strain typing, ANI > 99% + biology
- Species circumscription radius varies by genus (Parks 2018); 94-99% range
- GTDB taxonomy differs from NCBI for ~10% of species; document choice
- Type strain database matters; TYGS has built-in type-strain comparison
- 16S rRNA-based species ID is deprecated; do not rely on for new claims
- Always report AF alongside ANI; AF < 0.5 invalidates species call

## Related Skills

comparative-genomics/pangenome-analysis - ANI clustering precedes pangenome
comparative-genomics/ortholog-inference - Orthology consistency check
comparative-genomics/hgt-detection - HGT context for distant ANI
phylogenetics/species-trees - Phylogenomic species tree complements ANI
metagenomics/kraken-classification - Metagenomic classification
metagenomics/metaphlan-profiling - Strain-level profiling
genome-assembly/assembly-qc - QC before ANI
variant-calling/clinical-interpretation - Pathogen typing context
