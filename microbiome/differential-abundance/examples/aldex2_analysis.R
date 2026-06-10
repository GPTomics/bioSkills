# Reference: ALDEx2 1.34+, ANCOMBC 2.4+ | Verify API if version differs
# Consensus differential abundance on an amplicon feature table.
# Runs ALDEx2 (conservative Dirichlet-MC CLR) AND LinDA (fast CLR regression) on the
# same synthetic count table, then reports the intersection (high-confidence) and union
# (exploratory) of the two significant-taxa sets - because the hit list depends more on
# the tool than on the biology (Nearing 2022 Nat Commun 13:342), a single tool is not a
# defensible deliverable.
library(ALDEx2)
library(MicrobiomeStat)

set.seed(1)
n_per_group <- 20
n_taxa <- 60

# Synthetic compositional counts: a Dirichlet-Multinomial per sample. Ten taxa are shifted
# between groups (the truth); the rest are null. Counts (not proportions) are what ALDEx2
# and LinDA expect.
base <- rep(5, n_taxa)
shifted <- base
shifted[1:10] <- shifted[1:10] * 4
draw_sample <- function(alpha) {
    p <- as.numeric(rgamma(length(alpha), alpha, 1))
    p <- p / sum(p)
    rmultinom(1, size = 30000, prob = p)[, 1]
}
control <- sapply(1:n_per_group, function(i) draw_sample(base))
treated <- sapply(1:n_per_group, function(i) draw_sample(shifted))
counts <- cbind(control, treated)
rownames(counts) <- paste0('ASV', seq_len(n_taxa))
colnames(counts) <- paste0('S', seq_len(2 * n_per_group))
groups <- factor(rep(c('control', 'treated'), each = n_per_group))
meta <- data.frame(Group = groups, row.names = colnames(counts))

# prev_cut 0.10: a feature must appear in >= 10% of samples. A modeling knob, not housekeeping -
# it reshapes the BH denominator. Raising to 0.25 removes more tests (more power on survivors,
# fewer rare-but-real taxa). Declare it and confirm results are not knife-edge-sensitive.
prev_cut <- 0.10
keep <- rowSums(counts > 0) >= prev_cut * ncol(counts)
counts <- counts[keep, ]
cat('Taxa after prevalence filter:', nrow(counts), '\n')

# q_cut 0.05: BH false-discovery threshold across taxa. Uncorrected p is meaningless with
# dozens-thousands of features. effect_floor 1: ALDEx2 effect is a standardized median-ratio;
# |effect| > 1 is a strong ~2-SD signal (Gloor 2016). Gate on q AND effect, not p alone.
q_cut <- 0.05
effect_floor <- 1.0

# mc_samples 128: Monte-Carlo Dirichlet draws; ALDEx2 reports the EXPECTED p over draws, which
# is why it is conservative. Use 256+ for publication. denom 'all' = CLR against the geometric
# mean of every feature (the default reference frame).
aldex_out <- aldex(counts, as.character(groups), mc.samples = 128, test = 't',
                   effect = TRUE, denom = 'all')
# we.eBH = Welch expected BH-adjusted p (report this, NOT we.ep)
sig_aldex <- rownames(aldex_out)[aldex_out$we.eBH < q_cut & abs(aldex_out$effect) > effect_floor]
cat('ALDEx2 significant taxa:', length(sig_aldex), '\n')

# LinDA: CLR regression with mode-based bias correction; asymptotic FDR; no Monte-Carlo so it
# scales. A random effect in the formula (e.g. + (1|SubjectID)) would make it a mixed model for
# repeated measures. feature.dat.type 'count' lets LinDA add its own pseudocount before CLR.
linda_out <- linda(feature.dat = as.data.frame(counts), meta.dat = meta,
                   formula = '~ Group', feature.dat.type = 'count',
                   prev.filter = 0, alpha = q_cut)
# names(linda_out$output) follow the model-matrix coefficient columns; index by position to be
# robust to the exact factor-level suffix. reject = padj <= alpha after BH correction.
linda_res <- linda_out$output[[1]]
sig_linda <- rownames(linda_res)[linda_res$reject]
cat('LinDA significant taxa:', length(sig_linda), '\n')

# Consensus: intersect the significant SETS (never pool p-values across tools).
confident <- intersect(sig_aldex, sig_linda)
exploratory <- union(sig_aldex, sig_linda)
cat('\nConsensus (intersection, high-confidence):', length(confident), '\n')
cat('  ', paste(sort(confident), collapse = ', '), '\n')
cat('Union (exploratory):', length(exploratory), '\n')

consensus <- data.frame(
    taxon = sort(exploratory),
    in_aldex = sort(exploratory) %in% sig_aldex,
    in_linda = sort(exploratory) %in% sig_linda
)
consensus$n_tools <- consensus$in_aldex + consensus$in_linda
consensus <- consensus[order(-consensus$n_tools, consensus$taxon), ]
write.csv(consensus, 'da_consensus.csv', row.names = FALSE)
cat('\nSaved: da_consensus.csv (taxon / which tools / agreement count)\n')
