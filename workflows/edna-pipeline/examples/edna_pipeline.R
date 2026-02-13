# Reference: DADA2 1.30+, FastQC 0.12+, MultiQC 1.21+, cutadapt 4.4+, phyloseq 1.46+, vegan 2.6+ | Verify API if version differs
## eDNA metabarcoding pipeline (DADA2 path): FASTQ to community ecology.
## Marker: COI (Leray primers) - adjust truncLen and reference DB for other markers.

library(dada2)
library(decontam)
library(phyloseq)
library(iNEXT)
library(vegan)
library(indicspecies)

# --- Configuration ---
FASTQ_DIR <- 'trimmed/'
METADATA_FILE <- 'metadata.csv'
OUTPUT_PREFIX <- 'edna_results'

# truncLen: set from quality profiles; marker-dependent
# COI Leray ~313bp: c(250, 200); 12S MiFish ~170bp: c(140, 130)
# ITS: do NOT truncate (variable length); use maxEE only
TRUNC_LEN <- c(250, 200)
# maxEE c(2,2): standard; relax to c(5,5) for degraded eDNA (e.g., sediment cores)
MAX_EE <- c(2, 2)
# minOverlap 20: standard for most eDNA markers; increase for short amplicons
MIN_OVERLAP <- 20
# Reference database: marker-specific (BOLD/Midori2 for COI, UNITE for ITS, SILVA for 16S/18S)
TAX_DB <- 'reference_db.fa.gz'
SPECIES_DB <- 'species_db.fa.gz'
# minBoot 50: sensitive taxonomy assignment; 80 for conservative
MIN_BOOT <- 50
# decontam prevalence threshold 0.5: standard; lower to 0.1 for stringent
DECONTAM_THRESHOLD <- 0.5
# Tag-jump filter: 0.1% of max abundance per ASV
TAG_JUMP_THRESHOLD <- 0.001

# --- Step 1: Load and inspect reads ---
meta <- read.csv(METADATA_FILE)

fnFs <- sort(list.files(FASTQ_DIR, pattern = '_R1.fastq.gz', full.names = TRUE))
fnRs <- sort(list.files(FASTQ_DIR, pattern = '_R2.fastq.gz', full.names = TRUE))
sample_names <- gsub('_R1.fastq.gz', '', basename(fnFs))

message(sprintf('Loaded: %d samples', length(sample_names)))

# --- Step 2: Filter and trim ---
filtFs <- file.path('filtered', paste0(sample_names, '_F_filt.fastq.gz'))
filtRs <- file.path('filtered', paste0(sample_names, '_R_filt.fastq.gz'))
dir.create('filtered', showWarnings = FALSE)

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen = TRUNC_LEN, maxEE = MAX_EE,
                     minLen = 50, truncQ = 2, rm.phix = TRUE,
                     multithread = TRUE)
message(sprintf('Reads passing filter: %d / %d (%.1f%%)',
                sum(out[, 2]), sum(out[, 1]),
                sum(out[, 2]) / sum(out[, 1]) * 100))

# QC gate: reads per sample >1000
low_samples <- sample_names[out[, 2] < 1000]
if (length(low_samples) > 0) {
    message(sprintf('WARNING: %d samples with <1000 reads: %s',
                    length(low_samples), paste(low_samples, collapse = ', ')))
}

# --- Step 3: Learn error rates ---
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# --- Step 4: Denoise ---
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

# --- Step 5: Merge paired reads ---
merged <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap = MIN_OVERLAP)

seqtab <- makeSequenceTable(merged)
message(sprintf('ASVs before chimera removal: %d', ncol(seqtab)))

# --- Step 6: Remove chimeras ---
# method 'consensus': standard; 'pooled' for higher sensitivity with fewer samples
seqtab_nochim <- removeBimeraDenovo(seqtab, method = 'consensus', multithread = TRUE)
chimera_rate <- 1 - sum(seqtab_nochim) / sum(seqtab)
message(sprintf('Chimera rate: %.1f%% (%d ASVs remaining)',
                chimera_rate * 100, ncol(seqtab_nochim)))

# QC gate: chimera rate <20%
if (chimera_rate > 0.20) {
    message('WARNING: High chimera rate. Check primer removal completeness and PCR conditions.')
}

# --- Step 7: Assign taxonomy ---
taxa <- assignTaxonomy(seqtab_nochim, TAX_DB, multithread = TRUE, minBoot = MIN_BOOT)
taxa <- addSpecies(taxa, SPECIES_DB)

assigned_phylum <- sum(!is.na(taxa[, 'Phylum']))
assignment_rate <- assigned_phylum / nrow(taxa) * 100
message(sprintf('Taxonomy assignment rate (phylum): %.1f%% (%d / %d)',
                assignment_rate, assigned_phylum, nrow(taxa)))

# QC gate: assignment rate marker-specific (>90% for COI/BOLD, >60% for understudied)
if (assignment_rate < 60) {
    message('WARNING: Low taxonomy assignment. Check reference database completeness.')
}

# --- Step 8: Build phyloseq object ---
rownames(meta) <- meta$sample_id
ps <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows = FALSE),
               tax_table(taxa),
               sample_data(meta))
message(sprintf('Phyloseq: %d samples, %d ASVs', nsamples(ps), ntaxa(ps)))

# --- Step 9: Contamination filtering (decontam) ---
sample_data(ps)$is_neg <- sample_data(ps)$sample_type == 'negative_control'
contam <- isContaminant(ps, method = 'prevalence', neg = 'is_neg',
                        threshold = DECONTAM_THRESHOLD)
n_contam <- sum(contam$contaminant)
message(sprintf('Contaminant ASVs identified: %d', n_contam))

ps_clean <- prune_taxa(!contam$contaminant, ps)
ps_clean <- subset_samples(ps_clean, sample_type != 'negative_control')

# Tag-jumping removal: cross-contamination from index hopping
otu <- as(otu_table(ps_clean), 'matrix')
max_per_asv <- apply(otu, 2, max)
for (j in 1:ncol(otu)) {
    threshold <- max_per_asv[j] * TAG_JUMP_THRESHOLD
    otu[otu[, j] < threshold, j] <- 0
}
otu_table(ps_clean) <- otu_table(otu, taxa_are_rows = FALSE)

# Remove empty ASVs after filtering
ps_clean <- prune_taxa(taxa_sums(ps_clean) > 0, ps_clean)
message(sprintf('After decontamination: %d samples, %d ASVs', nsamples(ps_clean), ntaxa(ps_clean)))

# --- Step 10: Diversity analysis (iNEXT Hill numbers) ---
otu_for_inext <- as(otu_table(ps_clean), 'matrix')

# Hill numbers: q=0 (richness), q=1 (Shannon equivalent), q=2 (Simpson equivalent)
# endpoint: extrapolate to 2x observed sample size for fair comparison
inext_out <- iNEXT(as.list(as.data.frame(t(otu_for_inext))),
                   q = c(0, 1, 2), datatype = 'abundance',
                   endpoint = 2 * max(rowSums(otu_for_inext)))

completeness <- inext_out$DataInfo$SC
message(sprintf('Sample completeness range: %.1f%% - %.1f%%',
                min(completeness) * 100, max(completeness) * 100))

# QC gate: completeness >80%
if (min(completeness) < 0.80) {
    message('WARNING: Some samples have low completeness. Deeper sequencing recommended.')
}

# --- Step 11: Community comparison (vegan) ---
env_data <- as(sample_data(ps_clean), 'data.frame')

# Hellinger transformation: standard for species composition; reduces dominant species influence
otu_hell <- decostand(otu_for_inext, method = 'hellinger')

# DCA on untransformed data to determine gradient length
dca <- decorana(otu_for_inext)
gradient_length <- diff(range(scores(dca, display = 'sites', choices = 1)))
message(sprintf('DCA gradient length: %.2f SD', gradient_length))

# RDA: linear response (<=3 SD), uses Hellinger-transformed data
# CCA: unimodal response (>3 SD), uses raw abundances (chi-squared distance)
# Replace temperature + depth + season with actual column names from metadata
if (gradient_length <= 3) {
    ord <- rda(otu_hell ~ temperature + depth + season, data = env_data)
    method_name <- 'RDA'
} else {
    ord <- cca(otu_for_inext ~ temperature + depth + season, data = env_data)
    method_name <- 'CCA'
}

# permutations 999: standard; increase to 9999 for publication
anova_result <- anova.cca(ord, permutations = 999)
message(sprintf('%s significance: F = %.2f, p = %.4f',
                method_name, anova_result$F[1], anova_result$`Pr(>F)`[1]))

# Indicator species analysis
indval <- multipatt(otu_for_inext, env_data$site, control = how(nperm = 999))
sig_indicators <- indval$sign[indval$sign$p.value < 0.05, ]
message(sprintf('Significant indicator ASVs: %d', nrow(sig_indicators)))

# --- Save results ---
write.csv(as.data.frame(tax_table(ps_clean)), paste0(OUTPUT_PREFIX, '_taxonomy.csv'))
write.csv(as.data.frame(otu_table(ps_clean)), paste0(OUTPUT_PREFIX, '_otu_table.csv'))
write.csv(as.data.frame(sig_indicators), paste0(OUTPUT_PREFIX, '_indicators.csv'))
saveRDS(ps_clean, paste0(OUTPUT_PREFIX, '_phyloseq.rds'))
saveRDS(inext_out, paste0(OUTPUT_PREFIX, '_inext.rds'))

message(sprintf('\nPipeline complete: %d samples, %d ASVs, completeness %.0f-%.0f%%',
                nsamples(ps_clean), ntaxa(ps_clean),
                min(completeness) * 100, max(completeness) * 100))
