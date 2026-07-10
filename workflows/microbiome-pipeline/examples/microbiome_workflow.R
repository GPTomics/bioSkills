# Reference: DADA2 1.30+, phyloseq 1.46+, ALDEx2 1.34+ | Verify API if version differs
# End-to-end 16S amplicon workflow: demultiplexed FASTQ -> consensus differential abundance.
# Orchestration only: per-step decisions live in the six microbiome category skills.
# Stage order: primers (cutadapt, BEFORE truncation) -> per-RUN DADA2 -> mergeSequenceTables ->
# one chimera removal -> taxonomy -> phyloseq + placed tree -> declared-depth diversity (adonis2 +
# betadisper) -> consensus DA (ALDEx2 + ANCOM-BC2) on UNRAREFIED counts.
library(dada2)
library(phyloseq)
library(vegan)
library(ALDEx2)
library(ANCOMBC)
library(ggplot2)

raw_path <- 'raw_reads'
silva_train <- 'silva_nr99_v138.1_train_set.fa.gz'
silva_species <- 'silva_species_assignment_v138.1.fa.gz'
metadata_file <- 'sample_metadata.csv'
placed_tree_file <- 'sepp_tree.nwk'   # SEPP/Greengenes2 PLACED tree; NOT a de novo tree from short reads
output_dir <- 'microbiome_results'
fwd_primer <- 'GTGYCAGCMGCCGCGGTAA'   # 515F
rev_primer <- 'GGACTACNVGGGTWTCTAAT'  # 806R
trunc_len <- c(240, 160)   # merge budget: truncLen_F + truncLen_R >= amplicon_len + ~12 (V4 ~253 bp has slack)
max_ee <- c(2, 2)          # expected-errors filter; reads above this are discarded (Callahan 2016)
min_boot <- 50             # assignTaxonomy bootstrap floor; ranks below it are NA (RDP, reads <=250 nt)
rare_depth <- 10000        # rarefaction depth on the alpha-rarefaction plateau; NOT min(sample_sums)
prev_cut <- 0.10           # DA prevalence filter: a feature must appear in >=10% of samples (declared knob)
mc_samples <- 128          # ALDEx2 Monte-Carlo draws; 256+ for publication
effect_floor <- 1          # ALDEx2 |effect| floor: between-group diff exceeds within-condition dispersion (Gloor 2016), gate WITH q

dir.create(output_dir, showWarnings = FALSE)
dir.create('trimmed', showWarnings = FALSE)
dir.create('filtered', showWarnings = FALSE)

# === 1. REMOVE PRIMERS (cutadapt, BEFORE any quality/error step) ===
# Leftover primers corrupt the per-run error model and masquerade as chimeras. -g/-G are the 5'
# forward/reverse primers; --discard-untrimmed drops primerless pairs (a primerless read is suspect).
fnFs_raw <- sort(list.files(raw_path, pattern = '_R1_001.fastq.gz', full.names = TRUE))
fnRs_raw <- sort(list.files(raw_path, pattern = '_R2_001.fastq.gz', full.names = TRUE))
sample_names <- sapply(strsplit(basename(fnFs_raw), '_'), `[`, 1)
fnFs <- file.path('trimmed', basename(fnFs_raw))
fnRs <- file.path('trimmed', basename(fnRs_raw))
for (i in seq_along(fnFs_raw)) {
    system2('cutadapt', c('-g', fwd_primer, '-G', rev_primer, '--discard-untrimmed',
                          '-o', fnFs[i], '-p', fnRs[i], fnFs_raw[i], fnRs_raw[i]))
}

# === 2-3. PER-RUN DADA2 -> mergeSequenceTables -> ONE chimera removal ===
# The error model is fit PER sequencing run. A run_id column in the metadata declares which run each
# sample belongs to; learnErrors runs once per run, never pooled. Single-run studies have one group.
metadata <- read.csv(metadata_file, row.names = 1)
run_of <- if ('run_id' %in% colnames(metadata)) metadata[sample_names, 'run_id'] else rep('run1', length(sample_names))

getN <- function(x) sum(getUniques(x))

denoise_one_run <- function(run_samples) {
    filtFs <- file.path('filtered', paste0(run_samples, '_F_filt.fastq.gz'))
    filtRs <- file.path('filtered', paste0(run_samples, '_R_filt.fastq.gz'))
    names(filtFs) <- run_samples
    names(filtRs) <- run_samples
    idx <- match(run_samples, sample_names)
    # Capture the filterAndTrim return: it is the ONLY source of per-sample input/filtered counts,
    # and without it the read-tracking checkpoint below cannot be computed.
    out <- filterAndTrim(fnFs[idx], filtFs, fnRs[idx], filtRs, truncLen = trunc_len, maxEE = max_ee,
                         truncQ = 2, maxN = 0, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
    errF <- learnErrors(filtFs, multithread = TRUE)   # this run only
    errR <- learnErrors(filtRs, multithread = TRUE)
    dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
    dadaRs <- dada(filtRs, err = errR, multithread = TRUE)
    mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
    # With ONE sample in a run, dada() returns a bare dada object and mergePairs() a bare data.frame,
    # not lists -- sapply would then iterate their internal slots/columns. Handle that case separately.
    getN_all <- function(x) if (length(run_samples) == 1) getN(x) else sapply(x, getN)
    track <- cbind(out, getN_all(dadaFs), getN_all(dadaRs), getN_all(mergers))
    colnames(track) <- c('input', 'filtered', 'denoisedF', 'denoisedR', 'merged')
    rownames(track) <- run_samples
    list(seqtab = makeSequenceTable(mergers), track = track)
}

run_out <- lapply(unique(run_of), function(r) denoise_one_run(sample_names[run_of == r]))
run_tables <- lapply(run_out, `[[`, 'seqtab')
seqtab <- if (length(run_tables) > 1) do.call(mergeSequenceTables, run_tables) else run_tables[[1]]
seqtab_nochim <- removeBimeraDenovo(seqtab, method = 'consensus', multithread = TRUE, verbose = TRUE)

# Per-sample read tracking through filter -> denoise -> merge -> nonchim. An aggregate percentage
# hides a per-sample merge cliff: a sample can retain 90% overall while losing most reads at merging
# because truncLen left no overlap. Inspect the MERGED column per sample, not just the total.
track <- do.call(rbind, lapply(run_out, `[[`, 'track'))
track <- cbind(track, nonchim = rowSums(seqtab_nochim)[rownames(track)])
print(track)
# >20% loss at the merge step is the conventional cliff threshold: it almost always means truncLen
# left <12 nt overlap (amplicon too long for the read pair), not a biological signal.
merge_loss_pct <- 100 * (1 - track[, 'merged'] / track[, 'denoisedF'])
if (any(merge_loss_pct > 20, na.rm = TRUE)) {
    warning('Merge cliff in: ', paste(rownames(track)[merge_loss_pct > 20], collapse = ', '),
            ' -- re-check truncLen leaves >=12 nt overlap plus biological length variation')
}
cat('ASVs after chimera removal:', ncol(seqtab_nochim),
    ' reads retained:', round(100 * sum(seqtab_nochim) / sum(seqtab), 1), '%\n')

# === 4. TAXONOMY (region-matched reference; genus for 16S, not species) ===
taxa <- assignTaxonomy(seqtab_nochim, silva_train, minBoot = min_boot, tryRC = TRUE, multithread = TRUE)
if (file.exists(silva_species)) {
    taxa <- addSpecies(taxa, silva_species)   # species ONLY by exact amplicon match, else NA
}

# === 5. phyloseq + PLACED tree (SEPP/Greengenes2, NOT de novo from short reads) ===
ps <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows = FALSE), tax_table(taxa), sample_data(metadata))
if (file.exists(placed_tree_file)) {
    ps <- merge_phyloseq(ps, ape::read.tree(placed_tree_file))   # tree tips must match the ASV-sequence taxa_names
}
asv_ids <- setNames(paste0('ASV', seq(ntaxa(ps))), taxa_names(ps))   # rename tree tips and taxa in lockstep
if (!is.null(phy_tree(ps, errorIfNULL = FALSE))) {
    phy_tree(ps)$tip.label <- asv_ids[phy_tree(ps)$tip.label]
}
taxa_names(ps) <- unname(asv_ids)
ps <- subset_taxa(ps, is.na(Order) | Order != 'Chloroplast')      # universal 16S primers amplify host organelle rRNA
ps <- subset_taxa(ps, is.na(Family) | Family != 'Mitochondria')   # remove before diversity/DA or closure distorts both
saveRDS(ps, file.path(output_dir, 'phyloseq_object.rds'))

# === 6. DIVERSITY (rarefy ONLY here; declare depth; report dropped samples) ===
below_depth <- names(which(sample_sums(ps) < rare_depth))
cat('Sampling depth:', rare_depth, ' samples dropped below it:', length(below_depth),
    if (length(below_depth)) paste0('(', paste(below_depth, collapse = ', '), ')') else '', '\n')
ps_rare <- rarefy_even_depth(ps, sample.size = rare_depth, rngseed = 42, replace = FALSE)

alpha_div <- estimate_richness(ps_rare, measures = c('Observed', 'Shannon', 'InvSimpson'))
alpha_div$Shannon_eff <- exp(alpha_div$Shannon)   # effective species (Hill q=1), base-invariant
alpha_div <- cbind(alpha_div, sample_data(ps_rare))
write.csv(alpha_div, file.path(output_dir, 'alpha_diversity.csv'))

# Prefer weighted UniFrac when a tree is present; fall back to Bray-Curtis otherwise.
beta_dist <- if (!is.null(phy_tree(ps_rare, errorIfNULL = FALSE))) UniFrac(ps_rare, weighted = TRUE) else phyloseq::distance(ps_rare, method = 'bray')
meta_df <- data.frame(sample_data(ps_rare))
permanova <- adonis2(beta_dist ~ Group, data = meta_df, permutations = 999)
dispersion <- permutest(betadisper(beta_dist, meta_df$Group))   # MANDATORY: location vs dispersion
cat('PERMANOVA R2:', round(permanova$R2[1], 3), ' p:', permanova$`Pr(>F)`[1],
    ' betadisper p:', dispersion$tab$`Pr(>F)`[1], '\n')

pcoa <- ordinate(ps_rare, method = 'PCoA', distance = beta_dist)
p_beta <- plot_ordination(ps_rare, pcoa, color = 'Group') + stat_ellipse(level = 0.95) + theme_minimal() +
    labs(title = sprintf('PCoA (PERMANOVA R2=%.2f, p=%.3f; betadisper p=%.3f)',
                         permanova$R2[1], permanova$`Pr(>F)`[1], dispersion$tab$`Pr(>F)`[1]))
ggsave(file.path(output_dir, 'beta_diversity_pcoa.pdf'), p_beta, width = 7, height = 6)

# === 7. CONSENSUS DIFFERENTIAL ABUNDANCE (UNRAREFIED counts; >=2 CoDA tools) ===
ps_filt <- filter_taxa(ps, function(x) sum(x > 0) >= prev_cut * nsamples(ps), TRUE)
counts <- as.matrix(otu_table(ps_filt))
if (!taxa_are_rows(ps_filt)) {
    counts <- t(counts)   # ALDEx2 wants integer counts with taxa in ROWS
}
groups <- as.character(sample_data(ps_filt)$Group)

# Multi-run study: carry run as a batch covariate. ANCOM-BC2 -> fix_formula = 'run_id + Group';
# ALDEx2 cannot take a covariate via aldex() -- use aldex.clr() + aldex.glm() on a model matrix.
ax <- aldex(counts, groups, mc.samples = mc_samples, test = 't', effect = TRUE, denom = 'all')
sig_aldex <- rownames(ax)[ax$we.eBH < 0.05 & abs(ax$effect) > effect_floor]   # q AND effect, not p alone

ab <- ancombc2(data = ps_filt, fix_formula = 'Group', p_adj_method = 'BH',   # single-run; multi-run -> 'run_id + Group'
               prv_cut = prev_cut, group = 'Group', struc_zero = TRUE, pseudo_sens = TRUE)$res
diff_col <- grep('^diff_Group', colnames(ab), value = TRUE)[1]
ss_col <- grep('^passed_ss_Group', colnames(ab), value = TRUE)[1]
sig_ancombc <- ab$taxon[ab[[diff_col]] & ab[[ss_col]]]   # significant AND pseudo-count-robust

confident <- intersect(sig_aldex, sig_ancombc)   # high-confidence; union = exploratory
exploratory <- union(sig_aldex, sig_ancombc)
consensus <- data.frame(ASV = exploratory,
                        aldex = exploratory %in% sig_aldex,
                        ancombc = exploratory %in% sig_ancombc,
                        confident = exploratory %in% confident)
write.csv(consensus, file.path(output_dir, 'da_consensus.csv'), row.names = FALSE)
cat('DA hits - ALDEx2:', length(sig_aldex), ' ANCOM-BC2:', length(sig_ancombc),
    ' confident (intersection):', length(confident), '\n')

cat('Pipeline complete. Results in', output_dir, '\n')
