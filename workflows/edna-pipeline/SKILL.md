---
name: bio-workflows-edna-pipeline
description: End-to-end eDNA metabarcoding from raw amplicons to community ecology. Covers QC, primer removal, denoising with OBITools3 or DADA2, contamination filtering, taxonomy assignment, Hill number diversity, and constrained ordination. Use when processing environmental DNA samples for biodiversity assessment or ecological surveys.
tool_type: mixed
primary_tool: obitools3
workflow: true
depends_on:
  - ecological-genomics/edna-metabarcoding
  - ecological-genomics/biodiversity-metrics
  - ecological-genomics/community-ecology
  - read-qc/quality-reports
qc_checkpoints:
  - after_demux: "Reads per sample >1000; negative controls <100 reads"
  - after_denoising: "Chimera rate <20%; ASV/OTU count reasonable for marker and environment"
  - after_decontam: "No unexpected taxa remaining in negative controls; tag-jumping artifacts removed"
  - after_taxonomy: "Assignment rate marker-specific: >90% to phylum for COI/BOLD, >60% for understudied markers"
  - after_diversity: "Rarefaction curves approaching asymptote; sample completeness >80%"
---

# eDNA Metabarcoding Pipeline

Complete workflow from raw amplicon sequences to community ecology analysis, supporting
both OBITools3 and DADA2 processing paths.

## Pipeline Overview

```
Raw amplicon FASTQ (demultiplexed)
    |
    v
[1. QC] ------------------> FastQC / MultiQC quality assessment
    |
    v
[2. Primer Removal] ------> Cutadapt (remove forward + reverse primers)
    |                            |
    |                            +---> QC: reads per sample >1000
    |
    +--- Path A: OBITools3     +--- Path B: DADA2
    |       |                  |       |
    |       v                  |       v
    |   [3a. obi alignpairedend]   [3b. filterAndTrim]
    |       |                  |       |
    |       v                  |       v
    |   [4a. obi uniq]        |   [4b. learnErrors + dada]
    |       |                  |       |
    |       v                  |       v
    |   [5a. obi ecotag]      |   [5b. assignTaxonomy]
    |                          |
    +-------- Merge -----------+
                |
                v
[6. Contamination Filter] -> decontam / microDecon (negative control removal)
    |
    v
[7. Taxonomy Table] -------> Species x sample matrix
    |
    v
[8. Diversity Analysis] ---> iNEXT Hill numbers (q=0,1,2)
    |
    v
[9. Community Comparison] -> vegan CCA/RDA + indicspecies
    |
    v
Species table + diversity metrics + ordination plots
```

## Step 1: Quality Assessment

```bash
fastqc -t 8 -o fastqc_output/ raw_reads/*.fastq.gz
multiqc fastqc_output/ -o multiqc_report/
```

## Step 2: Primer Removal

```bash
# Adapter sequences are marker-specific; examples below for common eDNA markers
# --discard-untrimmed: remove reads without primers (likely off-target)
# --minimum-length 50: discard very short fragments after trimming

# COI (Leray primers mlCOIintF / jgHCO2198)
cutadapt -g GGWACWGGWTGAACWGTWTAYCCYCC -G TAIACYTCIGGRTGICCRAARAAYCA \
    --discard-untrimmed --minimum-length 50 -j 8 \
    -o trimmed/{sample}_R1.fastq.gz -p trimmed/{sample}_R2.fastq.gz \
    raw_reads/{sample}_R1.fastq.gz raw_reads/{sample}_R2.fastq.gz
```

Common primer sets by marker:

| Marker | Forward Primer | Reverse Primer | Target |
|--------|---------------|----------------|--------|
| COI | mlCOIintF | jgHCO2198 | Metazoan invertebrates |
| 12S MiFish | MiFish-U-F | MiFish-U-R | Fish |
| ITS2 | ITS3 | ITS4 | Fungi |
| rbcL | rbcLa-F | rbcLa-R | Plants |
| 18S V9 | 1389F | 1510R | Eukaryotes |

### QC Checkpoint: Demultiplexing

```bash
# Gate: reads per sample >1000; negative controls <100 reads
for f in trimmed/*_R1.fastq.gz; do
    sample=$(basename "$f" _R1.fastq.gz)
    count=$(zcat "$f" | awk 'END{print NR/4}')
    echo "$sample: $count reads"
done
```

## Step 3: Paired-End Merging and Denoising

### Path A: OBITools3

```bash
# Import paired FASTQ into OBITools3 DMS
obi import --fastq-input trimmed/reads_R1.fastq.gz reads/reads1
obi import --fastq-input trimmed/reads_R2.fastq.gz reads/reads2

# Paired-end alignment
obi alignpairedend -R reads/reads2 reads/reads1 reads/aligned

# Filter by alignment score
# score >= 50: removes poorly overlapping pairs
obi grep -p 'sequence["score"] >= 50' reads/aligned reads/filtered

# Filter by merged length (marker-dependent range)
# 100-500 bp: typical COI amplicon range
obi grep -p 'len(sequence) >= 100 and len(sequence) <= 500' \
    reads/filtered reads/length_filtered

# Dereplicate
obi uniq reads/length_filtered reads/dereplicated

# Remove singletons
# count >=2: removes sequencing errors; increase to 5-10 for noisy datasets
obi grep -p 'sequence["count"] >= 2' reads/dereplicated reads/denoised

# Denoise (remove PCR/sequencing errors)
# ratio 0.05: sequences <5% abundance of a 1-mismatch parent are merged
obi clean -s merged_sample -r 0.05 -H reads/denoised reads/cleaned
```

### Path B: DADA2 (R)

```r
library(dada2)

path <- 'trimmed/'
fnFs <- sort(list.files(path, pattern = '_R1.fastq.gz', full.names = TRUE))
fnRs <- sort(list.files(path, pattern = '_R2.fastq.gz', full.names = TRUE))
sample_names <- gsub('_R1.fastq.gz', '', basename(fnFs))

# Filter and trim
# truncLen: set based on quality profiles; marker-dependent
# maxEE c(2,2): max expected errors; standard for eDNA
# minLen 50: minimum after trimming
filtFs <- file.path('filtered', paste0(sample_names, '_F_filt.fastq.gz'))
filtRs <- file.path('filtered', paste0(sample_names, '_R_filt.fastq.gz'))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen = c(200, 180), maxEE = c(2, 2),
                     minLen = 50, truncQ = 2, rm.phix = TRUE,
                     multithread = TRUE)

# Learn error rates
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# Denoise
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

# Merge paired reads
# minOverlap 20: standard; increase if amplicon has short overlap region
merged <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap = 20)

# Build ASV table
seqtab <- makeSequenceTable(merged)

# Remove chimeras
# method 'consensus': standard; 'pooled' for higher sensitivity
seqtab_nochim <- removeBimeraDenovo(seqtab, method = 'consensus', multithread = TRUE)
chimera_rate <- 1 - sum(seqtab_nochim) / sum(seqtab)
message(sprintf('Chimera rate: %.1f%%', chimera_rate * 100))
```

### QC Checkpoint: Denoising

```r
# Gate 1: Chimera rate <20%
if (chimera_rate > 0.20) message('WARNING: High chimera rate. Check primer removal and PCR conditions.')

# Gate 2: ASV count reasonable for marker
n_asvs <- ncol(seqtab_nochim)
message(sprintf('ASVs after denoising: %d', n_asvs))
# Typical ranges: COI 500-5000, 12S 50-500, ITS 200-3000
```

## Step 4: Contamination Filtering

### R (decontam)

```r
library(decontam)
library(phyloseq)

ps <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows = FALSE),
               sample_data(meta))

# Identify negative controls
sample_data(ps)$is_neg <- sample_data(ps)$sample_type == 'negative_control'

# prevalence method: standard for eDNA; threshold 0.5 identifies contaminants
# present more in negative controls than real samples
contam <- isContaminant(ps, method = 'prevalence', neg = 'is_neg', threshold = 0.5)
message(sprintf('Contaminant ASVs: %d', sum(contam$contaminant)))

ps_clean <- prune_taxa(!contam$contaminant, ps)

# Remove negative control samples
ps_clean <- subset_samples(ps_clean, sample_type != 'negative_control')
```

### Tag-jumping removal

```r
# Tag-jumping: cross-contamination from index hopping during sequencing
# Remove ASVs with <0.1% of max abundance in a sample (likely tag-jump artifacts)
otu <- as(otu_table(ps_clean), 'matrix')
max_per_asv <- apply(otu, 2, max)
otu_filtered <- otu
for (j in 1:ncol(otu)) {
    # 0.1% of max: standard tag-jump threshold
    threshold <- max_per_asv[j] * 0.001
    otu_filtered[otu[, j] < threshold, j] <- 0
}
otu_table(ps_clean) <- otu_table(otu_filtered, taxa_are_rows = FALSE)
```

### QC Checkpoint: Decontamination

```r
# Gate: verify contaminant ASVs were removed from real samples
n_before <- ntaxa(ps)
n_after <- ntaxa(ps_clean)
message(sprintf('ASVs removed as contaminants: %d (%.1f%%)',
                n_before - n_after, (n_before - n_after) / n_before * 100))
if ((n_before - n_after) / n_before > 0.5) {
    message('WARNING: >50% ASVs flagged as contaminants. Review decontam threshold.')
}
```

## Step 5: Taxonomy Assignment

### OBITools3 (ecotag)

```bash
# ecotag assigns taxonomy using LCA algorithm against reference database
# Reference databases: EMBL, BOLD, MIDORI2, UNITE (marker-dependent)
obi ecotag -R reads/refdb --taxonomy reads/taxonomy reads/cleaned reads/assigned

# Filter by assignment quality (species-level for COI)
obi grep -p 'sequence["best_identity"] >= 0.97' reads/assigned reads/filtered_assigned

obi export --tab-output reads/filtered_assigned > taxonomy_results.tsv
```

### DADA2 (assignTaxonomy)

```r
# SILVA for 16S/18S, UNITE for ITS, custom for COI/12S
# Reference databases must be formatted for DADA2
# minBoot 50: minimum bootstrap confidence; 80 for more conservative assignments
taxa <- assignTaxonomy(seqtab_nochim, 'reference_db.fa.gz', multithread = TRUE, minBoot = 50)
taxa <- addSpecies(taxa, 'species_db.fa.gz')
```

### Marker-specific taxonomy databases

| Marker | Database | Typical Assignment Rate |
|--------|----------|------------------------|
| COI | BOLD / Midori2 | >90% to phylum, 60-80% to species |
| 12S | MitoFish / 12S-seqdb | >90% to family for fish |
| ITS | UNITE | >80% to genus for fungi |
| rbcL | GenBank / NCBI nt | >85% to family for plants |
| 18S | SILVA / PR2 | >90% to phylum |

### QC Checkpoint: Taxonomy

```r
# Gate: assignment rate should meet marker expectations
assigned <- !is.na(taxa[, 'Phylum'])
assignment_rate <- sum(assigned) / length(assigned) * 100
message(sprintf('Taxonomy assignment rate (phylum level): %.1f%%', assignment_rate))
if (assignment_rate < 60) message('WARNING: Low assignment rate. Check reference database completeness.')
```

## Step 6: Diversity Analysis

### R (iNEXT)

```r
library(iNEXT)

otu_matrix <- as(otu_table(ps_clean), 'matrix')

# Hill numbers: q=0 (richness), q=1 (Shannon diversity), q=2 (Simpson diversity)
# endpoint: extrapolation to 2x observed sample size
inext_out <- iNEXT(as.list(as.data.frame(t(otu_matrix))),
                   q = c(0, 1, 2), datatype = 'abundance',
                   endpoint = 2 * max(rowSums(otu_matrix)))

# Sample completeness: fraction of estimated diversity observed
completeness <- inext_out$DataInfo$SC
message(sprintf('Sample completeness range: %.1f%% - %.1f%%',
                min(completeness) * 100, max(completeness) * 100))
```

### QC Checkpoint: Diversity

```r
# Gate 1: rarefaction approaching asymptote
if (min(completeness) < 0.80) {
    message('WARNING: Some samples have low completeness (<80%). Deeper sequencing recommended.')
}

# Gate 2: reasonable richness estimates
q0_estimates <- inext_out$AsyEst[inext_out$AsyEst$Diversity == 'Species richness', ]
message(sprintf('Estimated richness range: %d - %d species',
                min(q0_estimates$Estimator), max(q0_estimates$Estimator)))
```

## Step 7: Community Comparison

### R (vegan + indicspecies)

```r
library(vegan)
library(indicspecies)

otu_matrix <- as(otu_table(ps_clean), 'matrix')
env_data <- as(sample_data(ps_clean), 'data.frame')

# Hellinger transformation: standard for community composition data
# Reduces influence of dominant species
otu_hell <- decostand(otu_matrix, method = 'hellinger')

# DCA on untransformed data to determine gradient length
dca <- decorana(otu_matrix)
gradient_length <- diff(range(scores(dca, display = 'sites', choices = 1)))
message(sprintf('DCA gradient length: %.2f SD', gradient_length))

# RDA: linear response (<=3 SD), uses Hellinger-transformed data
# CCA: unimodal response (>3 SD), uses raw abundances (chi-squared distance)
if (gradient_length <= 3) {
    ord <- rda(otu_hell ~ temperature + depth + season, data = env_data)
    method_name <- 'RDA'
} else {
    ord <- cca(otu_matrix ~ temperature + depth + season, data = env_data)
    method_name <- 'CCA'
}

# Permutation test for significance
# permutations 999: standard; increase to 9999 for publication
anova_result <- anova.cca(ord, permutations = 999)
message(sprintf('%s significance: p = %.4f', method_name, anova_result$`Pr(>F)`[1]))

# Indicator species analysis
# Groups defined by environmental category (e.g., site, season, habitat)
indval <- multipatt(otu_matrix, env_data$site, control = how(nperm = 999))
summary(indval, alpha = 0.05)
```

## Parameter Recommendations

| Step | Parameter | Recommendation |
|------|-----------|----------------|
| Cutadapt | --discard-untrimmed | Always use; removes off-target reads |
| Cutadapt | --minimum-length | 50 (general); adjust per expected amplicon size |
| DADA2 | truncLen | Set from quality profiles; marker-dependent |
| DADA2 | maxEE | c(2,2) standard; c(5,5) for degraded eDNA |
| DADA2 | minOverlap | 20 (standard); increase for short overlaps |
| OBITools3 | --min-count | 2 (removes singletons); 5-10 for noisy datasets |
| decontam | threshold | 0.5 (prevalence method); 0.1 for stringent |
| Tag-jumping | threshold | 0.1% of max abundance per ASV |
| Taxonomy | minBoot | 50 (sensitive); 80 (conservative) |
| ecotag | --min-identity | 0.97 (COI species); 0.95 (genus); marker-dependent |
| iNEXT | endpoint | 2x max observed sample size |
| vegan | permutations | 999 (standard); 9999 (publication) |

## Troubleshooting

| Issue | Likely Cause | Solution |
|-------|--------------|----------|
| Few reads after primer removal | Wrong primer sequences or orientation | Verify primer sequences; try --revcomp |
| High chimera rate (>20%) | Excessive PCR cycles or low-quality input | Reduce PCR cycles; improve DNA extraction |
| Many unassigned ASVs | Incomplete reference database | Use marker-specific database; lower minBoot |
| Contamination in negatives | Tag-jumping or lab contamination | Apply tag-jump filter; review extraction protocol |
| Low sample completeness | Insufficient sequencing depth | Increase sequencing; pool fewer samples |
| Ordination axes not significant | Weak environmental gradients | Add more environmental variables; check sample size |
| Unexpected taxa (e.g., human) | Sample contamination | Filter known contaminants; review field protocols |
| Very few ASVs | Over-aggressive filtering | Relax truncLen, maxEE, or min-count thresholds |

## Related Skills

- ecological-genomics/edna-metabarcoding - Detailed eDNA processing
- ecological-genomics/biodiversity-metrics - Diversity analysis details
- ecological-genomics/community-ecology - Ordination and indicator species
- read-qc/quality-reports - Raw read quality assessment
- microbiome/amplicon-processing - 16S clinical alternative
