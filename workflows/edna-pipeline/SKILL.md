---
name: bio-workflows-edna-pipeline
description: End-to-end eDNA metabarcoding from raw amplicons to community ecology. Covers QC, primer removal (mandatory before DADA2 filterAndTrim), denoising with OBITools3 v3 (obi stats plural; DMS-based) or DADA2 ASVs (Callahan 2017), decontam combined method as screening-not-classifier (Davis 2018), tag-jumping with NovaSeq 10x MiSeq caveat (Schnell 2015), Hill-number effective species counts with coverage-based rarefaction (Jost 2006; Chao & Jost 2012; doubling rule), beta-diversity decomposition with MANDATORY PERMANOVA + PERMDISP pair (Anderson & Walsh 2013), constrained ordination, and the read-counts-not-abundance critique (Lamb 2019). Use when processing eDNA samples for biodiversity assessment, deciding ASV vs OTU, configuring OBITools3 v3, interpreting decontam screening, or reporting community comparisons with the dispersion confound check.
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
  - after_denoising: "Chimera rate <20% (>30% indicates library-prep issues); ASV/OTU count reasonable for marker"
  - after_decontam: "decontam combined method at threshold 0.1 (0.05 for low-biomass); biological-plausibility review of each flagged ASV; tag-jumping rate quantified and filtered (~0.001-0.005 MiSeq; ~0.005-0.01 NovaSeq per Schnell 2015)"
  - after_taxonomy: "Assignment rate marker-specific: 50-85% unassigned is typical per Wangensteen 2018; report gap honestly"
  - after_diversity: "Hill numbers reported as effective species counts (not raw Shannon); coverage-based rarefaction at C=0.95; extrapolation bounded at 2x reference (doubling rule); sample completeness >80%"
  - after_ordination: "PERMANOVA + PERMDISP reported TOGETHER (Anderson & Walsh 2013); if betadisper significant, location conclusion is not supported"
---

## Version Compatibility

Reference examples tested with: DADA2 1.30+, FastQC 0.12+, MultiQC 1.21+, cutadapt 4.4+, phyloseq 1.46+, vegan 2.6+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# eDNA Metabarcoding Pipeline

**"Process my eDNA samples from raw reads to community ecology"** → Orchestrate primer removal, denoising (OBITools3 or DADA2), contamination filtering, taxonomy assignment, Hill number diversity estimation, and constrained ordination for species-environment analysis.

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

# Combined method (Davis 2018): uses BOTH negative controls AND DNA concentration
# threshold=0.1 default; 0.05 for low-biomass samples
# CRITICAL: output is SCREENING, not classification; biological-plausibility check required
if ('dna_concentration' %in% sample_variables(ps)) {
    contam <- isContaminant(ps, method = 'combined', neg = 'is_neg',
                            conc = 'dna_concentration', threshold = 0.1)
} else {
    # Fall back to prevalence-only if no qPCR/Qubit DNA-concentration data
    contam <- isContaminant(ps, method = 'prevalence', neg = 'is_neg',
                            threshold = 0.1)
}
message(sprintf('Flagged candidate contaminant ASVs: %d', sum(contam$contaminant)))
message('Manual review required: verify biological plausibility before deletion')

ps_clean <- prune_taxa(!contam$contaminant, ps)

# Remove negative control samples
ps_clean <- subset_samples(ps_clean, sample_type != 'negative_control')
```

### Tag-jumping removal (Schnell 2015; NovaSeq caveat)

```r
# Tag-jumping: cross-contamination from index hopping during library prep / sequencing
# Schnell 2015 Mol Ecol Resour 15:1289-1303 documented 0.1-2% per read pair
# NovaSeq patterned flow cells have ~10x higher rates than MiSeq
# Use platform-appropriate threshold:
#   MiSeq: 0.001-0.005 (0.1-0.5% of ASV total)
#   NovaSeq: 0.005-0.01 (0.5-1% of ASV total)
# Quantify residual rate first from per-ASV cross-sample appearance
otu <- as(otu_table(ps_clean), 'matrix')
max_per_asv <- apply(otu, 2, max)
otu_filtered <- otu
# Set threshold by platform; default below is for MiSeq
tag_jump_frac <- 0.001
for (j in 1:ncol(otu)) {
    threshold <- max_per_asv[j] * tag_jump_frac
    otu_filtered[otu[, j] < threshold, j] <- 0
}
otu_table(ps_clean) <- otu_table(otu_filtered, taxa_are_rows = FALSE)
# Modern alternative: metabaR::tagjumpslayer(metabarlist_obj, threshold = 0.03)
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

# MANDATORY companion: PERMANOVA + PERMDISP (Anderson & Walsh 2013 Ecol Monogr 83:557-574)
# adonis2 tests centroid difference; betadisper tests dispersion homogeneity
# If betadisper is also significant, PERMANOVA significance is dispersion-confounded
bray_dist <- vegdist(otu_matrix, method = 'bray')
permanova <- adonis2(bray_dist ~ site, data = env_data,
                     by = 'margin', permutations = 999)
disp <- betadisper(bray_dist, env_data$site)
disp_test <- permutest(disp, permutations = 999)
message(sprintf('PERMANOVA p = %.4f; PERMDISP p = %.4f',
                permanova[['Pr(>F)']][1], disp_test$tab[['Pr(>F)']][1]))
if (permanova[['Pr(>F)']][1] < 0.05 && disp_test$tab[['Pr(>F)']][1] < 0.05) {
    message('WARNING: Both PERMANOVA and PERMDISP significant; location-vs-dispersion confounded')
}

# Indicator species analysis with group-size equalization (NOT basic IndVal)
# func='IndVal.g' corrects for unbalanced group sizes (De Caceres & Legendre 2009)
indval <- multipatt(otu_matrix, env_data$site, func = 'IndVal.g',
                    control = how(nperm = 999))
summary(indval, alpha = 0.05)
```

## Parameter Recommendations

| Step | Parameter | Recommendation |
|------|-----------|----------------|
| Cutadapt | --discard-untrimmed | Always use; removes off-target reads. MANDATORY before DADA2 filterAndTrim |
| Cutadapt | --minimum-length | 50 (general); adjust per expected amplicon size |
| DADA2 | truncLen | Set from quality profiles; marker-dependent (typical 2x250 COI: c(220,180)) |
| DADA2 | maxEE | c(2,2) standard; c(5,5) for degraded eDNA |
| DADA2 | minOverlap | 20 (standard); increase for short overlaps |
| DADA2 | chimera method | 'consensus' standard (conservative); 'pooled' more aggressive |
| OBITools3 v3 | command syntax | `obi stats` (plural; was `obistat` in v1); `.tar.gz` taxonomy |
| OBITools3 | --min-count | 2 (removes singletons); 5-10 for noisy datasets |
| decontam | method | 'combined' if concentration AND controls; 'prevalence' fallback |
| decontam | threshold | 0.1 default; 0.05 for low-biomass samples |
| Tag-jumping MiSeq | threshold | 0.001-0.005 fraction of ASV total |
| Tag-jumping NovaSeq | threshold | 0.005-0.01 (~10x MiSeq per Schnell 2015) |
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
