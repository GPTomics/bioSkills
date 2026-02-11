# eDNA Metabarcoding Usage Guide

## Overview

Environmental DNA (eDNA) metabarcoding identifies species from trace DNA shed into water, soil, or air. This skill covers the full workflow from raw amplicon reads to species occurrence tables: primer removal, quality filtering, denoising (OBITools3 or DADA2), taxonomy assignment against curated reference databases (MIDORI2, MitoFish, BOLD), contamination filtering with decontam, and occupancy modeling for detection probability correction.

## Prerequisites

- OBITools3 (`pip install OBITools3`)
- cutadapt (`pip install cutadapt`)
- R with DADA2 (`BiocManager::install('dada2')`)
- decontam (`BiocManager::install('decontam')`)
- microDecon (`remotes::install_github('donaldtmcknight/microDecon')`)
- occumb (`install.packages('occumb')`) and JAGS for occupancy modeling
- Reference database formatted for the target marker (MIDORI2, MitoFish, BOLD)

## Quick Start

Tell your AI agent what you want to do:

- "Process my eDNA water samples to identify fish species using 12S MiFish primers"
- "Run the OBITools3 pipeline on my COI metabarcoding data"
- "Denoise my eDNA amplicon reads with DADA2 and assign taxonomy against MIDORI2"
- "Filter contamination from my eDNA dataset using extraction blanks"
- "Run occupancy modeling on my replicated eDNA samples to correct for detection probability"

## Example Prompts

### Basic Processing

> "I have paired-end FASTQ files from eDNA water samples amplified with MiFish 12S primers. Remove primers with cutadapt, denoise with DADA2, and assign taxonomy against MitoFish."

> "Run the full OBITools3 pipeline on my COI metabarcoding data: align paired ends, dereplicate, denoise, and assign taxonomy."

### Contamination Control

> "I have 48 eDNA samples and 6 extraction blanks. Use decontam to identify and remove contaminant ASVs using both frequency and prevalence methods."

> "My eDNA dataset has tag-jumping artifacts between samples. Use microDecon to clean up the cross-contamination."

### Occupancy Modeling

> "I collected 3 eDNA replicates per site across 20 sites. Fit an occupancy model with occumb to estimate true species occurrence probabilities."

### Database Preparation

> "Download and format the MIDORI2 COI database for DADA2 taxonomy assignment."

## What the Agent Will Do

1. Remove primer sequences with cutadapt using linked adapter mode for the target barcode region
2. Process reads through either OBITools3 (dereplication + denoising) or DADA2 (ASV inference)
3. Filter sequences by length and quality appropriate for the target marker
4. Remove chimeric sequences and PCR/sequencing artifacts
5. Assign taxonomy against the appropriate reference database with confidence thresholds
6. Filter contamination using negative controls (decontam) and tag-jumping correction (microDecon)
7. Optionally fit occupancy models to estimate detection-corrected species occurrence
8. Export a species-by-site occurrence or abundance table

## Tips

- Always include extraction blanks and PCR negative controls in sequencing runs for decontam filtering
- COI has a wider barcode gap than 12S, making species-level assignment more reliable for invertebrates
- MiFish 12S is preferred for fish eDNA due to shorter amplicon and universal primer performance
- DADA2 produces exact ASVs (single-nucleotide resolution) while OBITools3 clusters at a similarity threshold
- For rare species detection, aim for at least 50,000 reads per sample after filtering
- Biological replicates (separate water collections) are more informative than PCR replicates for occupancy modeling
- Always check chimera rates: above 20% suggests issues with PCR conditions or primer dimers
- The MIDORI2 database is updated regularly; always use the latest version for best taxonomy coverage

## Related Skills

- biodiversity-metrics - Diversity analysis from species occurrence tables
- community-ecology - Environmental gradient analysis of communities
- microbiome/amplicon-processing - 16S clinical microbiome alternative
- read-qc/quality-reports - Upstream read quality assessment
