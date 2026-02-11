# eDNA Metabarcoding Pipeline - Usage Guide

## Overview
Complete environmental DNA metabarcoding workflow from raw amplicon sequences to biodiversity assessment and community ecology. Supports two processing paths: OBITools3 (CLI-based, optimized for eDNA) and DADA2 (R-based, ASV resolution). Includes contamination filtering with decontam, Hill number diversity analysis with iNEXT, and constrained ordination with vegan for community comparison. Handles all common eDNA markers (COI, 12S, ITS, rbcL, 18S).

## Prerequisites
```bash
# OBITools3 (CLI path)
pip install obitools3
# or conda install -c bioconda obitools3

# General CLI tools
pip install cutadapt
# or conda install -c bioconda cutadapt fastqc multiqc

# R packages (DADA2 path + downstream analysis)
install.packages('BiocManager')
BiocManager::install(c('dada2', 'decontam', 'phyloseq'))
install.packages(c('iNEXT', 'vegan', 'indicspecies'))

# Reference databases (download once, marker-specific):
# COI: BOLD reference library or Midori2
# 12S: MitoFish or 12S-seqdb
# ITS: UNITE (https://unite.ut.ee/)
# 16S/18S: SILVA (https://www.arb-silva.de/)
```

**Input data:**
- Demultiplexed paired-end FASTQ files (one pair per sample)
- Sample metadata with sample type (field sample vs. negative control)
- Environmental variables for community analysis (e.g., temperature, depth, site)
- Primer sequences for the target marker

## Quick Start
Tell your AI agent what you want to do:
- "I have eDNA amplicon data from river water samples - process it from raw reads to species lists"
- "Run the full eDNA metabarcoding pipeline on my COI amplicon data with OBITools3"
- "Process my 12S MiFish eDNA data through DADA2 and calculate fish diversity"
- "Compare eDNA communities across sampling sites with constrained ordination"
- "Filter contamination from my eDNA dataset using negative controls"

## Example Prompts

### Full Pipeline
> "I have paired-end COI amplicon data from 50 water samples plus 5 negative controls. Process everything from raw FASTQ to species diversity estimates."

> "Run the complete eDNA pipeline: trim primers, denoise with DADA2, remove contaminants, assign taxonomy against BOLD, and calculate Hill number diversity."

### Specific Steps
> "I already have an ASV table from DADA2. Run contamination filtering with decontam using my negative controls, then assign taxonomy."

> "Compare fish communities detected by 12S eDNA across upstream and downstream sites using constrained ordination."

### Marker-Specific
> "Process my fungal ITS2 eDNA data using DADA2 with the UNITE database for taxonomy assignment."

> "I have eDNA samples from a marine survey using the MiFish 12S primers. Run the full pipeline optimized for fish detection."

## What the Agent Will Do
1. Run FastQC and MultiQC on raw reads to assess quality
2. Remove primer sequences with Cutadapt (marker-specific primers)
3. Validate: check reads per sample >1000, negative controls <100 reads
4. Merge paired ends and denoise (OBITools3 or DADA2 depending on preference)
5. Remove chimeric sequences
6. Validate: check chimera rate <20%, ASV count reasonable for marker
7. Filter contamination using negative controls (decontam prevalence method)
8. Remove tag-jumping artifacts
9. Validate: negative controls clean after filtering
10. Assign taxonomy using marker-appropriate reference database
11. Validate: assignment rate meets marker expectations
12. Calculate Hill number diversity (q=0, 1, 2) with iNEXT
13. Validate: rarefaction approaching asymptote, sample completeness >80%
14. Run constrained ordination (RDA or CCA) and indicator species analysis
15. Export species table, diversity metrics, and ordination plots

## Tips
- Always include negative controls (extraction blanks and PCR blanks) for contamination filtering
- OBITools3 is optimized for eDNA workflows and handles tag-jumping natively; DADA2 provides ASV-level resolution
- Primer removal is critical: untrimmed primers cause artificial diversity inflation
- Tag-jumping (index hopping) is a major concern for eDNA on Illumina platforms; always apply the 0.1% threshold filter
- Taxonomy assignment quality depends heavily on reference database completeness; COI/BOLD gives best species-level results for animals
- Hill numbers provide a unified diversity framework: q=0 (richness), q=1 (Shannon equivalent), q=2 (Simpson equivalent)
- Use DCA gradient length to choose between RDA (linear, <3 SD) and CCA (unimodal, >3 SD)
- For degraded eDNA (e.g., sediment cores), relax DADA2 maxEE to c(5,5) and reduce truncLen
- Increase sequencing depth for samples with low completeness rather than lowering quality thresholds

## Related Skills
- ecological-genomics/edna-metabarcoding - Detailed eDNA processing
- ecological-genomics/biodiversity-metrics - Diversity analysis details
- ecological-genomics/community-ecology - Ordination and indicator species
- read-qc/quality-reports - Raw read quality assessment
- microbiome/amplicon-processing - 16S clinical alternative
