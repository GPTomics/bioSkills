# Reference: limorhyde 1.0+, limma 3.50+ | Verify API if version differs
# LimoRhyde + limma differential rhythmicity on simulated two-condition data.
# Fits condition*(time_cos+time_sin) and tests the condition:time INTERACTION
# (differential RHYTHMICITY) separately from the condition MAIN effect (differential
# EXPRESSION). An explicit interaction model borrows strength across conditions and
# avoids the detect-then-Venn anti-pattern that overestimates rhythm reprogramming.

suppressPackageStartupMessages({library(limorhyde); library(limma)})

set.seed(1)

period <- 24                     # circadian period in hours; time is measured (ZT), not inferred pseudotime
zt <- seq(0, 44, by = 4)         # 12 timepoints spanning ~2 full 24h cycles: the minimum for a rhythm test
n_rep <- 2                       # DR estimates an INTERACTION, so it needs more replication than detection
conditions <- c('WT', 'KO')

meta <- expand.grid(zt = zt, rep = seq_len(n_rep), condition = conditions)
meta$condition <- factor(meta$condition, levels = c('WT', 'KO'))   # WT is the reference group
meta$sample <- sprintf('%s_ZT%02d_r%d', meta$condition, meta$zt, meta$rep)

# Generator: per-condition amplitude and phase encode the four differential-rhythmicity classes.
sim_gene <- function(base, amp_wt, amp_ko, phase_wt, phase_ko, noise = 0.3) {
  amp <- ifelse(meta$condition == 'WT', amp_wt, amp_ko)
  phase <- ifelse(meta$condition == 'WT', phase_wt, phase_ko)
  base + amp * cos(2 * pi * (meta$zt - phase) / period) + rnorm(nrow(meta), 0, noise)
}

genes <- list(
  rhythmic_both = sim_gene(8, 1.0, 1.0, 6, 6),    # same rhythm both conditions -> NOT differentially rhythmic
  loss_in_KO    = sim_gene(8, 1.0, 0.0, 6, 6),    # loss of rhythm in KO
  gain_in_KO    = sim_gene(8, 0.0, 1.0, 6, 6),    # gain of rhythm in KO
  phase_shift   = sim_gene(8, 1.0, 1.0, 6, 14),   # ~8h phase change in KO
  amp_reduced   = sim_gene(8, 1.2, 0.3, 6, 6),    # amplitude change in KO
  flat_both     = sim_gene(8, 0.0, 0.0, 6, 6),    # arrhythmic in both
  DE_only       = sim_gene(8, 0.0, 0.0, 6, 6) + ifelse(meta$condition == 'KO', 3, 0))  # mean shift, no rhythm change

expr <- do.call(rbind, genes)
colnames(expr) <- meta$sample

# limorhyde() decomposes measured time into a cosinor basis; prefix 'time_' names them time_cos, time_sin.
meta <- cbind(meta, limorhyde(meta$zt, 'time_', period = period))

# Differential RHYTHMICITY = condition:time interaction; differential EXPRESSION = condition main effect.
design <- model.matrix(~ condition * (time_cos + time_sin), data = meta)
fit <- eBayes(lmFit(expr, design))

dr_cols <- grep('conditionKO:time_', colnames(design), value = TRUE)  # the two interaction coefficients
de_col <- 'conditionKO'                                               # condition main effect, adjusting for time

# Moderated F over BOTH interaction terms ranks differential rhythmicity; topTable returns BH-adjusted P.
dr <- topTable(fit, coef = dr_cols, number = Inf, sort.by = 'F')
de <- topTable(fit, coef = de_col, number = Inf, sort.by = 'p')

cat('Differential RHYTHMICITY (condition:time interaction), ranked:\n')
print(round(dr[, c('F', 'P.Value', 'adj.P.Val')], 4))
cat('\nDifferential EXPRESSION (condition main effect), ranked:\n')
print(round(de[, c('logFC', 'P.Value', 'adj.P.Val')], 4))

# Expected: loss/gain/phase/amp genes top the DR list; DE_only tops DE but is NOT differentially rhythmic;
# rhythmic_both and flat_both are neither. Classify DR direction from per-condition amplitude/phase estimates.
