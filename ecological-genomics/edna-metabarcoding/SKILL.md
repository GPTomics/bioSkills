---
name: bio-ecological-genomics-edna-metabarcoding
description: Processes environmental DNA metabarcoding data from raw amplicon reads to species occurrence tables using OBITools3, DADA2, and taxonomic assignment against BOLD, MIDORI2, or MitoFish databases. Handles COI, 12S, rbcL, and ITS barcode regions with primer removal, denoising, chimera detection, and contamination filtering via decontam. Includes occupancy modeling (occumb) for detection probability correction. Use when analyzing eDNA from water, soil, or bulk samples for biodiversity monitoring. Not for 16S human microbiome (see microbiome/amplicon-processing).
tool_type: mixed
primary_tool: obitools3
---

# eDNA Metabarcoding

Processes environmental DNA amplicon reads into species occurrence tables with taxonomy assignment, contamination filtering, and occupancy modeling.

## Primer Removal with Cutadapt

Linked adapter trimming removes primer pairs while discarding reads lacking primers:

```bash
# COI primers (mlCOIintF / jgHCO2198)
cutadapt -g 'GGWACWGGWTGAACWGTWTAYCCYCC;min_overlap=20' \
         -G 'TAIACYTCIGGRTGICCRAARAAYCA;min_overlap=20' \
         --discard-untrimmed --pair-filter=any \
         -o trimmed_R1.fastq.gz -p trimmed_R2.fastq.gz \
         raw_R1.fastq.gz raw_R2.fastq.gz

# 12S MiFish-U primers
cutadapt -g 'GTCGGTAAAACTCGTGCCAGC;min_overlap=18' \
         -G 'CATAGTGGGGTATCTAATCCCAGTTTG;min_overlap=18' \
         --discard-untrimmed --pair-filter=any \
         -o trimmed_R1.fastq.gz -p trimmed_R2.fastq.gz \
         raw_R1.fastq.gz raw_R2.fastq.gz

# ITS primers (ITS1F / ITS2)
cutadapt -g 'CTTGGTCATTTAGAGGAAGTAA;min_overlap=18' \
         -G 'GCTGCGTTCTTCATCGATGC;min_overlap=18' \
         --discard-untrimmed --pair-filter=any \
         -o trimmed_R1.fastq.gz -p trimmed_R2.fastq.gz \
         raw_R1.fastq.gz raw_R2.fastq.gz
```

## OBITools3 Pipeline

Full pipeline from paired-end reads to taxonomy table:

```bash
# Import paired FASTQ into OBITools3 DMS
obi import --fastq-input raw_R1.fastq.gz EDNA/reads1
obi import --fastq-input raw_R2.fastq.gz EDNA/reads2

# Paired-end alignment
obi alignpairedend -R EDNA/reads2 EDNA/reads1 EDNA/aligned

# Filter by alignment score and length
# score >= 50: removes poorly overlapping pairs
# length 100-500: typical COI amplicon range
obi grep -p 'sequence["score"] >= 50' EDNA/aligned EDNA/filtered
obi grep -p 'len(sequence) >= 100 and len(sequence) <= 500' \
    EDNA/filtered EDNA/length_filtered

# Demultiplex if needed (with ngsfilter file)
obi ngsfilter -t ngsfilter.txt -u EDNA/unassigned \
    EDNA/length_filtered EDNA/demux

# Dereplicate (obi uniq creates merged_sample attribute automatically)
obi uniq EDNA/demux EDNA/derep

# Remove singletons (count >= 2 removes PCR/sequencing errors)
obi grep -p 'sequence["count"] >= 2' EDNA/derep EDNA/no_singletons

# Denoise (remove PCR/sequencing errors)
# obi clean: denoises by merging low-abundance variants into parent sequences
# -s merged_sample: per-sample denoising (critical for multi-sample datasets)
# -r 0.05: ratio threshold; -H: keep only head sequences (discard variants)
obi clean -s merged_sample -r 0.05 -H EDNA/no_singletons EDNA/denoised

# Taxonomy assignment against reference database
obi ecotag -R EDNA/refdb --taxonomy EDNA/taxonomy EDNA/denoised EDNA/assigned

# Export to tab-separated format
obi export --tab-output EDNA/assigned > species_table.tsv
```

## DADA2 Pipeline for eDNA

DADA2 provides an alternative ASV-based approach:

```r
library(dada2)

# After cutadapt primer removal, get file paths
filt_fwd <- file.path('filtered', basename(fwd_reads))
filt_rev <- file.path('filtered', basename(rev_reads))

# Filter and trim
# maxEE=2: expected errors threshold balancing sensitivity/specificity
# truncLen: set based on quality profile inspection
out <- filterAndTrim(fwd_reads, filt_fwd, rev_reads, filt_rev,
                     maxN = 0, maxEE = c(2, 2), truncQ = 2,
                     minLen = 100, rm.phix = TRUE, multithread = TRUE)

# Learn error rates
err_fwd <- learnErrors(filt_fwd, multithread = TRUE)
err_rev <- learnErrors(filt_rev, multithread = TRUE)

# Denoise
dada_fwd <- dada(filt_fwd, err = err_fwd, multithread = TRUE)
dada_rev <- dada(filt_rev, err = err_rev, multithread = TRUE)

# Merge paired ends
# minOverlap=12: default; increase for short amplicons
merged <- mergePairs(dada_fwd, filt_fwd, dada_rev, filt_rev, minOverlap = 12)

# Build sequence table and remove chimeras
seqtab <- makeSequenceTable(merged)
# Chimera rate <20% is typical; >30% suggests library prep issues
seqtab_nochim <- removeBimeraDenovo(seqtab, method = 'consensus', multithread = TRUE)

# Taxonomy assignment against MIDORI2 or custom database
# minBoot=80: bootstrap confidence for genus-level; 50 for family
taxa <- assignTaxonomy(seqtab_nochim, 'MIDORI2_DADA2_COI.fasta.gz', minBoot = 80)
```

## Reference Databases

| Database | Markers | Format | Source |
|----------|---------|--------|--------|
| MIDORI2 | COI, 12S, 16S, 18S | DADA2, BLAST, OBITools | midori2.info |
| MitoFish | 12S (fish) | DADA2, BLAST | mitofish.aori.u-tokyo.ac.jp |
| BOLD | COI | FASTA, API | boldsystems.org |
| SILVA | 18S | DADA2, mothur | silva-project.net |
| UNITE | ITS (fungi) | DADA2, BLAST | unite.ut.ee |

Formatting a custom FASTA for DADA2 taxonomy assignment:

```bash
# DADA2 format: >Kingdom;Phylum;Class;Order;Family;Genus;Species
# Download MIDORI2 DADA2-formatted reference:
wget https://reference-midori.info/download/Databases/MIDORI2_DADA2/COI/MIDORI2_LONGEST_NUC_GB259_CO1_DADA2.fasta.gz
```

## Contamination Filtering

### decontam (frequency/prevalence method)

```r
library(decontam)

# Frequency method: uses DNA concentration from extraction blanks
contam_freq <- isContaminant(seqtab_nochim, conc = dna_concentration, method = 'frequency')

# Prevalence method: compares true samples vs negative controls
contam_prev <- isContaminant(seqtab_nochim, neg = is_negative_control, method = 'prevalence',
                             threshold = 0.5)
# threshold=0.5: ASV must be more prevalent in negatives than samples to flag

# Combined method
contam_both <- isContaminant(seqtab_nochim, conc = dna_concentration,
                             neg = is_negative_control, method = 'both')

seqtab_clean <- seqtab_nochim[, !contam_both$contaminant]
```

### microDecon for tag-jumping artifacts

```r
library(microDecon)

# Remove tag-jumping contamination between samples
# numb.blanks: number of blank/negative control columns at start of table
# numb.ind: vector of sample counts per group
decon_result <- decon(otu_table, numb.blanks = 3, numb.ind = c(10, 10, 10))
cleaned_table <- decon_result$decon.table
```

## Occupancy Modeling with occumb

Corrects for imperfect detection in metabarcoding replicates:

```r
library(occumb)

# Requires replicated eDNA samples (biological or PCR replicates)
# y: 3D array [species x sites x replicates] of read counts
data <- occumbData(y = count_array, spec_cov = species_covariates,
                   site_cov = site_covariates)

# Fit multi-species occupancy model (requires JAGS)
# n.chains=4, n.iter=10000: standard MCMC settings
fit <- occumb(data = data, n.chains = 4, n.iter = 10000, n.thin = 5, n.burn = 2500)

# Extract corrected occupancy estimates
summary(fit)
```

## Key Thresholds

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Min read depth per sample | 1,000 | Below this, rare species detection drops significantly |
| Singleton removal | count >= 2 | Singletons often represent sequencing errors |
| Chimera rate | <20% | Higher rates suggest library prep problems |
| COI species identity | >= 97% | Standard COI barcode gap threshold |
| COI genus identity | >= 95% | Conservative genus-level assignment |
| 12S species identity | >= 98% | 12S has less divergence than COI |
| decontam threshold | 0.5 | Balanced false positive/negative rate |
| DADA2 minBoot | 80 | Genus-level confidence; use 50 for family |

## Related Skills

- biodiversity-metrics - Diversity analysis from species occurrence tables
- community-ecology - Environmental gradient analysis of communities
- microbiome/amplicon-processing - 16S clinical microbiome alternative
- read-qc/quality-reports - Upstream read quality assessment
